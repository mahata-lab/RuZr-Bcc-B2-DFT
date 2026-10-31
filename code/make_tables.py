"""LaTeX tables for the manuscript (Table 2) and Supplementary Material (Tables S2-S5). Run after ruzr_analysis.py."""
import json, re, numpy as np, pandas as pd, os
import importlib.util, sys
spec = importlib.util.spec_from_file_location('ra', 'code/ruzr_analysis.py'); ra = importlib.util.module_from_spec(spec); spec.loader.exec_module(ra)
os.makedirs('results/latex', exist_ok=True)
D = pd.read_csv('results/configurations.csv'); C = pd.read_csv('results/compositions.csv', index_col=0)
DS = {e['id']: e for e in json.load(open('data/dataset.json'))}
famname = {'base':'base','X_for_Zr':'Z','X_for_Ru':'R','Ru_lean_X_rich':'R$^{*}$'}
order = {'base':0,'X_for_Zr':1,'X_for_Ru':2,'Ru_lean_X_rich':3}
stname = {'near-bcc':'near-BCC','omega-collapsed':'$\\omega$-collapsed','strongly distorted':'distorted','B off-site':'B off-site'}
def tex(lbl): return '$' + re.sub(r'([A-Z][a-z]?)(\d+)', r'\\mathrm{\1}_{\2}', lbl) + '$'
D['o'] = D.family.map(order); D = D.sort_values(['o','VEC','composition','configuration'])
def cfg(r): return tex(r.composition) + (f' ({r.configuration})' if (D.composition == r.composition).sum() > 1 else '')
# S2 convergence
with open('results/latex/table_convergence.tex','w') as f:
    for _, r in D.iterrows():
        f.write(f"{cfg(r)} & {famname[r.family]} & {int(r.ionic_steps)} & {r.max_force:.3f} & {r.pressure_kB:+.1f} & {r.sigma0_GPa:+.2f} \\\\\n")
# S3 descriptors
st = {'near-bcc':'n','omega-collapsed':'$\\omega$','strongly distorted':'d','B off-site':'B$^\\dagger$'}
with open('results/latex/table_configs.tex','w') as f:
    for _, r in D.iterrows():
        nef = '--' if np.isnan(r.NEF_tot) else f'{r.NEF_tot:.2f}'
        f.write(f"{cfg(r)} & {famname[r.family]} & {r.VEC:.3f} & {r.a_bcc:.3f} & {r.C11:.1f} & {r.C12:.1f} & {r.C44:.1f} & {r.Cp:.1f} & {r.B:.1f} & {r.G:.1f} & {r.G_over_B:.3f} & {r.nu:.3f} & {r.Cauchy:.1f} & {r.u_rms:.2f} & {r.omega_amp:.2f} & {nef} & {st[r.state]} \\\\\n")
# S4 tensors
with open('results/latex/table_tensors.tex','w') as f:
    for _, r in D.iterrows():
        Cm, _ = ra.elastic_tensor(DS[r.id]); f.write(f"\\multicolumn{{7}}{{l}}{{{cfg(r)}}}\\\\\n")
        for i in range(6): f.write(' & ' + ' & '.join(f'{Cm[i,j]:.1f}' for j in range(6)) + ' \\\\\n')
        f.write('\\addlinespace\n')
# S5 omega statistics
with open('results/latex/table_omega.tex','w') as f:
    for _, r in D.sort_values('omega_frac', ascending=False).iterrows():
        f.write(f"{cfg(r)} & {100*r.omega_frac:.1f} & {r.omega_random_percentile:.1f} & {r.kappa:+.2f} \\\\\n")
# Table 2 (main text)
assess = {('X_for_Zr','Re'):'$+$',('X_for_Zr','Mo'):'$+$ (weaker)',('X_for_Zr','Nb'):'$\\circ$ (lower misfit)',('X_for_Zr','Hf'):'$\\circ$',('X_for_Zr','Al'):'$-$ ($\\omega$)',
 ('X_for_Ru','Al'):"$-$ ($C'$)",('X_for_Ru','Sc'):"$-$ ($C'$)",('X_for_Ru','V'):"$-$ ($C'$)",('X_for_Ru','Y'):"$-$ ($C'$, misfit)",('X_for_Ru','Ta'):"$-$ ($C'$)",('X_for_Ru','Ti'):'$-$',('X_for_Ru','Re'):'$\\circ$',
 ('X_for_Ru','Mo'):'$-$ ($\\omega$)',('X_for_Ru','Nb'):'$-$ ($\\omega$)',('X_for_Ru','Hf'):'$-$ ($\\omega$)',('Ru_lean_X_rich','Ta'):"$-$ ($C'$)",('Ru_lean_X_rich','Ti'):'$-$ (distortion)',('Ru_lean_X_rich','Sc'):'$-$ (distortion)',('Ru_lean_X_rich','Y'):'$-$ (distortion)'}
groups = [('X_for_Zr',['Re','Hf','Al'],'\\emph{X for Zr at 7.4 at.\\% Ru (\\cmp{Zr48Ru4X2})}'),
          ('X_for_Ru',['Al','Sc','V','Y','Ti','Hf'],'\\emph{X for Ru, 3.7 at.\\% Ru (\\cmp{Zr50Ru2X2})}'),
          ('Ru_lean_X_rich',['Ta','Ti','Y'],'\\emph{Ru-lean, X-rich (\\cmp{Zr48Ru2X4})}')]
lines = []
for fam, xs, head in groups:
    lines.append('\\multicolumn{8}{l}{' + head + '}\\\\')
    for x in xs:
        r = C[(C.family == fam) & (C.X == x)].iloc[0]
        dc = f'{r.d_Cp:+.1f}'; dc = f'(${dc}$)' if r.state != 'near-bcc' else f'${dc}$'
        lines.append(f"{x} & {stname[r.state]} & {dc} & ${r.d_u_rms:+.3f}$ & ${r.d_omega_amp:+.2f}$ & ${r.d_G_over_B:+.3f}$ & ${r.d_misfit_pp:+.2f}$ & {assess[(fam,x)]} \\\\")
    lines.append('\\midrule')
open('results/latex/table_design.tex','w').write('\n'.join(lines[:-1]) + '\n')
print('tables written to results/latex/')

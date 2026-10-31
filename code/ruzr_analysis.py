"""
Analysis of the Zr-Ru(-X) BCC SQS calculations (54-atom cells).

Reads  data/dataset.json  (one entry per SQS configuration) and data/dos/*.npz, and writes results/*.csv:
  configurations.csv     elastic constants, moduli, lattice parameter, static displacement, omega projection,
                         N(E_F), classification (near-BCC / omega-collapsed / strongly distorted / B off-site)
  compositions.csv       composition means and changes relative to Zr50Ru4
  trajectories.csv       energy, RMS displacement and omega amplitude along each relaxation
  zr_local_environment.csv   displacement of every Zr site vs number of Ru in its first two shells
  fits.txt               linear fits and correlation coefficients quoted in the paper
Usage:  python code/ruzr_analysis.py          (run from the repository root)
"""
import json, numpy as np, pandas as pd
from scipy import stats

VALENCE = {'Zr':4,'Ru':8,'B':3,'Al':3,'Sc':3,'Y':3,'Ti':4,'Hf':4,'V':5,'Nb':5,'Ta':5,'Mo':6,'Re':7}
A_B2_EXP = 3.2505          # measured B2 ZrRu lattice parameter (Polash et al., J. Alloys Compd. 2026)
VARIANTS = [np.array(v) for v in ([1,1,1],[-1,1,1],[1,-1,1],[1,1,-1])]

# ---------------------------------------------------------------- elastic
def voigt_from_vasp(s):                       # VASP order xx yy zz xy yz zx -> Voigt xx yy zz yz xz xy
    xx,yy,zz,xy,yz,zx = s; return np.array([xx,yy,zz,yz,zx,xy])
def strain_voigt(L0, L):
    F = np.linalg.solve(L0, L).T; e = 0.5*(F+F.T) - np.eye(3)
    return np.array([e[0,0],e[1,1],e[2,2],2*e[1,2],2*e[0,2],2*e[0,1]])
def elastic_tensor(entry):
    L0 = np.array(entry['relaxed_lattice'])
    eps = np.array([strain_voigt(L0, np.array(s['lattice'])) for s in entry['strains']])
    sig = np.array([-0.1*voigt_from_vasp(s['stress_kB_vasp_order_xx_yy_zz_xy_yz_zx']) for s in entry['strains']])  # GPa, tension positive
    Amat = np.hstack([np.ones((len(eps),1)), eps]); coef = np.linalg.lstsq(Amat, sig, rcond=None)[0]
    C = coef[1:].T; return 0.5*(C+C.T), coef[0]
def moduli(C):
    S = np.linalg.inv(C)
    BV = (C[0,0]+C[1,1]+C[2,2]+2*(C[0,1]+C[0,2]+C[1,2]))/9
    GV = (C[0,0]+C[1,1]+C[2,2]-(C[0,1]+C[0,2]+C[1,2])+3*(C[3,3]+C[4,4]+C[5,5]))/15
    BR = 1/(S[0,0]+S[1,1]+S[2,2]+2*(S[0,1]+S[0,2]+S[1,2]))
    GR = 15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[0,2]+S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))
    B, G = (BV+BR)/2, (GV+GR)/2
    C11 = np.mean(np.diag(C)[:3]); C12 = np.mean([C[0,1],C[0,2],C[1,2]]); C44 = np.mean(np.diag(C)[3:])
    return dict(C11=C11, C12=C12, C44=C44, Cp=(C11-C12)/2, Cauchy=C12-C44, B=B, G=G, E=9*B*G/(3*B+G),
                nu=(3*B-2*G)/(2*(3*B+G)), G_over_B=G/B, Zener=C44/((C11-C12)/2), AU=5*GV/GR+BV/BR-6,
                Cmin_eig=np.linalg.eigvalsh(C).min())

# ---------------------------------------------------------------- displacements / omega
def displacement_field(L, X):
    n = np.round(np.asarray(X)*6).astype(int); u = ((np.asarray(X)*6-n)/6) @ np.asarray(L); return n, u-u.mean(0)
def omega_projection(n, u):
    tot = np.sum(u**2); best = None
    for v in VARIANTS:
        ph = 2*np.pi/3*(n@v); e = v/np.linalg.norm(v); ul = u@e
        c2, s2 = np.sum(np.cos(ph)**2), np.sum(np.sin(ph)**2)
        f = ((ul@np.cos(ph))**2/c2 + (ul@np.sin(ph))**2/s2)/tot
        c, s = ul@np.cos(ph)/c2, ul@np.sin(ph)/s2; phi0 = np.arctan2(-c, s)
        if best is None or f > best[0]: best = (f, np.cos(3*phi0))
    return best                                 # (f_omega, kappa)
def random_baseline(n, samples=20000, seed=42):
    rng = np.random.default_rng(seed); out = []
    for _ in range(samples):
        u = rng.normal(size=(len(n),3)); u -= u.mean(0); out.append(omega_projection(n,u)[0])
    return np.array(out)

# ---------------------------------------------------------------- DOS
def smooth(y, dE, s=0.10):
    k = np.arange(-4*s, 4*s+dE, dE); g = np.exp(-k**2/2/s**2); return np.convolve(y, g/g.sum(), 'same')
def dos_descriptors(path, species, counts):
    d = np.load(path); E = d['E']; dE = E[1]-E[0]; i0 = np.argmin(abs(E)); N = sum(counts)
    out = dict(NEF_tot=smooth(d['tot'],dE)[i0]/N)
    for el, c in zip(species, counts):
        out[f'NEF_d_{el}'] = smooth(d[el+'_t2g']+d[el+'_eg'], dE)[i0]/c
    return out

# ---------------------------------------------------------------- main
def classify(r):
    if r.composition.endswith('B2'): return 'B off-site'
    if r.omega_frac > 0.4: return 'omega-collapsed'
    if r.u_rms > 0.38: return 'strongly distorted'
    return 'near-bcc'

def main():
    ds = json.load(open('data/dataset.json')); rows=[]; traj=[]; local=[]
    base_n = None
    for e in ds:
        sp, cnt = e['species'], e['counts']; N = sum(cnt); comp = dict(zip(sp,cnt))
        L = np.array(e['relaxed_lattice']); V = abs(np.linalg.det(L))
        C, sig0 = elastic_tensor(e); m = moduli(C)
        n, u = displacement_field(L, e['relaxed_positions_direct']); f, kappa = omega_projection(n, u)
        lab = np.array(sum([[s]*c for s,c in zip(sp,cnt)],[])); um = np.linalg.norm(u,axis=1)
        X = [s for s in sp if s not in ('Zr','Ru')]
        r = dict(id=e['id'], composition=e['composition'], family=e['family'], configuration=e['configuration'],
                 X=X[0] if X else '', VEC=sum(VALENCE[s]*c for s,c in comp.items())/N, x_Ru=comp['Ru']/N,
                 a_bcc=(2*V/N)**(1/3), **m, sigma0_GPa=float(np.mean(sig0[:3])),
                 u_rms=float(np.sqrt(np.mean(um**2))), omega_frac=f, omega_amp=float(np.sqrt(f)*np.sqrt(np.mean(um**2))), kappa=kappa,
                 urms_Zr=float(np.sqrt(np.mean(um[lab=='Zr']**2))), urms_Ru=float(np.sqrt(np.mean(um[lab=='Ru']**2))),
                 ionic_steps=e['relax']['ionic_steps'], max_force=e['relax']['max_residual_force_eV_per_A'],
                 pressure_kB=e['relax']['pressure_kB'])
        if e['dos_file']: r.update(dos_descriptors(e['dos_file'], sp, cnt))
        rows.append(r)
        if base_n is None and e['composition']=='Zr50Ru4': base_n = n
        # local environment of Zr sites
        for i in np.where(lab=='Zr')[0]:
            dd = ((n-n[i]+3)%6)-3; d2 = (dd**2).sum(1)
            local.append(dict(id=e['id'], nRu=int(((lab=='Ru')&((d2==3)|(d2==4))).sum()), nRu_first=int(((lab=='Ru')&(d2==3)).sum()), u=um[i]))
        if e['relaxation_trajectory']:
            for k, st in enumerate(e['relaxation_trajectory']):
                nn, uu = displacement_field(st['lattice'], st['positions']); ff,_ = omega_projection(nn, uu)
                ur = float(np.sqrt(np.mean(np.sum(uu**2,1)))); traj.append(dict(id=e['id'], step=k, energy_eV=st['energy_eV'], u_rms=ur, omega_amp=float(np.sqrt(ff))*ur))
    D = pd.DataFrame(rows); D['state'] = D.apply(classify, axis=1)
    rb = random_baseline(base_n); D['omega_random_percentile'] = [(rb < x).mean()*100 for x in D.omega_frac]
    D.to_csv('results/configurations.csv', index=False)
    pd.DataFrame(traj).to_csv('results/trajectories.csv', index=False)
    Lc = pd.DataFrame(local); Lc.to_csv('results/zr_local_environment.csv', index=False)
    # composition means and changes relative to Zr50Ru4
    num = ['VEC','a_bcc','C11','C12','C44','Cp','B','G','G_over_B','nu','Cauchy','u_rms','omega_amp','NEF_tot']
    Cm = D.groupby('composition').agg({**{k:'mean' for k in num}, 'family':'first', 'X':'first', 'state':'first', 'id':'count'}).rename(columns={'id':'n_configurations'})
    b = Cm.loc['Zr50Ru4']
    for k in ['Cp','G_over_B','u_rms','omega_amp']: Cm['d_'+k] = Cm[k]-b[k]
    Cm['misfit_pct'] = 200*(Cm.a_bcc-A_B2_EXP)/(Cm.a_bcc+A_B2_EXP); Cm['d_misfit_pp'] = Cm.misfit_pct-Cm.loc['Zr50Ru4','misfit_pct']
    Cm.to_csv('results/compositions.csv')
    # fits quoted in the paper
    nb = D[D.state=='near-bcc']; W = D.dropna(subset=['NEF_tot'])
    with open('results/fits.txt','w') as f:
        def lin(x,y,data,tag):
            s = stats.linregress(data[x], data[y]); f.write(f'{tag}: slope {s.slope:.4f} +- {s.stderr:.4f}, r {s.rvalue:.3f}, p {s.pvalue:.2e}, n {len(data)}\n')
        lin('VEC','Cp',nb,"C' vs VEC (near-BCC configurations)")
        nbc = nb.groupby('composition')[['VEC','Cp']].mean(); lin('VEC','Cp',nbc,"C' vs VEC (near-BCC composition means)")
        lin('VEC','G_over_B',nb,'G/B vs VEC (near-BCC)')
        lin('VEC','NEF_tot',W[W.state=='near-bcc'],'N(EF) vs VEC (near-BCC)')
        lin('NEF_tot','Cp',W[W.state=='near-bcc'],"C' vs N(EF) (near-BCC)")
        lin('NEF_tot','G_over_B',W,'G/B vs N(EF) (all with DOS)')
        lin('G_over_B','Cauchy',D,'Cauchy vs G/B (all)')
        g = Lc.groupby('id').apply(lambda d: d[d.nRu==0].u.mean()-d[d.nRu>=1].u.mean())
        f.write(f'Zr local environment: Ru-free minus Ru-coordinated displacement > 0 in {(g>0).sum()} of {len(g)} configurations, '
                f'sign test p {stats.binomtest(int((g>0).sum()),len(g)).pvalue:.1e}, mean {g.mean():.3f} A\n')
        f.write(f'random-field f_omega: mean {rb.mean():.4f}, p99 {np.percentile(rb,99):.4f}, p99.9 {np.percentile(rb,99.9):.4f}\n')
    print(D[['id','state','Cp','G_over_B','u_rms','omega_amp','NEF_tot']].round(3).to_string()); print(open('results/fits.txt').read())

if __name__ == '__main__':
    main()

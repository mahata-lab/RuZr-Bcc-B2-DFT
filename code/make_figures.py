"""Figures 1-6 and Supplementary Figures S1-S2. Run from the repository root after code/ruzr_analysis.py."""
import pandas as pd, numpy as np, json, matplotlib, os
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.linewidth':0.8,'xtick.direction':'in','ytick.direction':'in',
 'xtick.top':True,'ytick.right':True,'legend.frameon':False,'axes.labelsize':9,'savefig.dpi':600,'lines.linewidth':1.6})
OUT='figures/'; os.makedirs(OUT,exist_ok=True)
FAMKEY={'base':'BASE','X_for_Zr':'Z','X_for_Ru':'R','Ru_lean_X_rich':'ZR'}
U=pd.read_csv('results/configurations.csv').rename(columns={'composition':'label','omega_frac':'omega_frac_max'})
U['fam_true']=U.family.map(FAMKEY); U['X']=U.X.fillna('')
DS={e['id']:e for e in json.load(open('data/dataset.json'))}
FAM={'BASE':('#0b0b0b','*',10,'Zr$_{50}$Ru$_4$ (base)'),'R':('#2a78d6','o',6,'Zr$_{50}$Ru$_2$X$_2$ (X for Ru)'),
     'Z':('#eb6834','s',6,'Zr$_{48}$Ru$_4$X$_2$ (X for Zr)'),'ZR':('#1baf7a','^',7,'Zr$_{48}$Ru$_2$X$_4$')}
def el(r): return 'base' if r.fam_true=='BASE' else r.X
def short(r):
    if r.fam_true=='BASE': return 'base'
    return {'R':r.X,'Z':r.X,'ZR':'4'+r.X}[r.fam_true]
U['short']=U.apply(short,axis=1)
def scatter(ax,x,y,data=U,label=True,states=True,fs=6.5,offs=None,skip=()):
    for f,(c,m,s,lab) in FAM.items():
        d=data[data.fam_true==f]
        for st,face in [('near-bcc',True),('B off-site',True),('omega-collapsed',False),('strongly distorted',False)]:
            e=d[d.state==st]
            if not len(e): continue
            ax.plot(e[x],e[y],m,ms=s,mfc=c if face else 'white',mec=c,mew=1.3 if not face else 0.6,ls='none',
                    alpha=0.55 if st=='B off-site' else 1)
        if label:
            for _,r in d.iterrows():
                if (r.fam_true,r.short) in skip: continue
                dx,dy=(offs or {}).get(r.short,(3,2))
                ax.annotate(r.short,(r[x],r[y]),xytext=(dx,dy),textcoords='offset points',fontsize=fs,color='#3a3a38')
def famlegend(ax,loc='best',extra=True):
    from matplotlib.lines import Line2D
    h=[Line2D([],[],marker=m,ls='none',mfc=c,mec=c,ms=s*0.8,label=l) for f,(c,m,s,l) in FAM.items()]
    if extra:
        h+=[Line2D([],[],marker='o',ls='none',mfc='white',mec='#52514e',mew=1.3,ms=5,label='ω-collapsed / strongly distorted')]
    ax.legend(handles=h,loc=loc,fontsize=6.5,handletextpad=0.3,borderaxespad=0.3)

# ---------- Fig 1: composition space + lattice parameter ----------
fig,ax=plt.subplots(1,2,figsize=(7.1,2.9),gridspec_kw={'width_ratios':[1,1.35]})
a=ax[0]
a.axhspan(8.0,9.7,color='#e8e6df',lw=0); a.text(0.3,8.85,'measured BCC matrix\n(Polash et al., 2026)',fontsize=6.3,va='center',color='#52514e')
pts={'BASE':(0,4/54*100),'R':(2/54*100,2/54*100),'Z':(2/54*100,4/54*100),'ZR':(4/54*100,2/54*100)}
for f,(xx,yy) in pts.items():
    c,m,s,l=FAM[f]; a.plot(xx,yy,m,ms=s+2,mfc=c,mec=c)
lbl={'BASE':'Zr$_{50}$Ru$_4$','R':'Zr$_{50}$Ru$_2$X$_2$\n(11 X)','Z':'Zr$_{48}$Ru$_4$X$_2$\n(Al,B,Hf,Mo,Nb,Re)','ZR':'Zr$_{48}$Ru$_2$X$_4$\n(Sc,Ta,Ti,Y)'}
off={'BASE':(0.25,0.25),'R':(0.25,-0.9),'Z':(0.25,0.25),'ZR':(0.25,-0.9)}
for f,(xx,yy) in pts.items(): a.text(xx+off[f][0],yy+off[f][1],lbl[f],fontsize=6.5)
a.annotate('',xy=(2/54*100,2/54*100+0.15),xytext=(0.1,4/54*100-0.15),arrowprops=dict(arrowstyle='->',color='#9a9890',lw=0.8))
a.annotate('',xy=(2/54*100-0.1,4/54*100),xytext=(0.25,4/54*100),arrowprops=dict(arrowstyle='->',color='#9a9890',lw=0.8))
a.set_xlim(-0.5,10.5); a.set_ylim(2.5,10.5); a.set_xlabel('X content (at.%)'); a.set_ylabel('Ru content in BCC (at.%)')
a.text(-0.2,1.04,'(a)',transform=a.transAxes,fontweight='bold')
b=ax[1]
C=U.groupby('label').agg(a_bcc=('a_bcc','mean'),fam=('fam_true','first'),short=('short','first'),state=('state','first')).reset_index()
C=C.sort_values('a_bcc')
cols=[FAM[f][0] for f in C.fam]
b.bar(range(len(C)),C.a_bcc,color=cols,width=0.75,edgecolor='white',linewidth=0.6)
b.axhline(3.5446,color='#0b0b0b',ls='--',lw=0.9); b.text(len(C)-0.5,3.5452,'exp. BCC (Zr$_{92}$Ru$_8$) 3.5446 Å',ha='right',va='bottom',fontsize=6.3)
b.set_xticks(range(len(C))); b.set_xticklabels(C.short,rotation=90,fontsize=6.5); b.set_ylim(3.50,3.61); b.set_ylabel('$a_\\mathrm{BCC}$ (Å)')
aB2=3.2505
sec=b.secondary_yaxis('right',functions=(lambda A:200*(A-aB2)/(A+aB2),lambda d:aB2*(200+d)/(200-d)))
sec.set_ylabel('BCC/B2 misfit δ (%)'); sec.tick_params(direction='in')
famlegend(b,loc='upper left',extra=False)
b.text(-0.12,1.04,'(b)',transform=b.transAxes,fontweight='bold')
fig.tight_layout(); fig.savefig(OUT+'Fig1_composition_lattice.png'); fig.savefig(OUT+'Fig1_composition_lattice.pdf'); plt.close()

# ---------- Fig 2: relaxation & omega ----------
T=pd.read_csv('results/trajectories.csv').rename(columns={'id':'tag','energy_eV':'E'})
paths={'base (config. 1)':'Zr50Ru4_c1','base (config. 2)':'Zr50Ru4_c2','Zr$_{48}$Ru$_4$Re$_2$':'Zr48Ru4Re2_c1','Zr$_{50}$Ru$_2$Hf$_2$':'Zr50Ru2Hf2_c1'}
pc=['#0b0b0b','#7a7974','#eb6834','#2a78d6']
fig,ax=plt.subplots(2,2,figsize=(7.1,5.2))
for (k,t),c in zip(paths.items(),pc):
    g=T[T.tag==t].sort_values('step'); E=(g.E-g.E.iloc[-1])*1000/54
    ax[0,0].plot(g.step,E,color=c,label=k); ax[0,1].plot(g.u_rms,g.omega_amp,color=c,label=k)
    ax[0,1].plot(g.u_rms.iloc[-1],g.omega_amp.iloc[-1],'o',color=c,ms=4)
ax[0,0].set_yscale('symlog',linthresh=1); ax[0,0].set_xlabel('ionic step'); ax[0,0].set_ylabel('$E-E_\\mathrm{final}$ (meV/atom)'); ax[0,0].legend(fontsize=6.5)
ax[0,1].plot([0,0.42],[0,0.42],':',color='#9a9890',lw=0.8); ax[0,1].text(0.27,0.29,'pure ω displacement',rotation=35,fontsize=6,color='#7a7974')
ax[0,1].axhline(np.sqrt(2/3)*np.sqrt(3)/12*3.55,color='#9a9890',lw=0.8,ls='--'); ax[0,1].text(0.005,0.405,'ideal ω (full collapse)',fontsize=6,color='#7a7974')
ax[0,1].set_xlabel('RMS static displacement $\\bar u$ (Å)'); ax[0,1].set_ylabel('ω-mode amplitude $u_\\omega$ (Å)'); ax[0,1].set_xlim(0,0.42); ax[0,1].set_ylim(0,0.45)
# plane resolved for Hf@Ru best variant and base cfg1
def plane(tag):
    r=DS[tag]; L=np.array(r['relaxed_lattice']); X=np.array(r['relaxed_positions_direct']); n=np.round(X*6).astype(int); u=((X*6-n)/6)@L; u-=u.mean(0)
    best=None
    for v in ([1,1,1],[-1,1,1],[1,-1,1],[1,1,-1]):
        v=np.array(v); e=v/np.linalg.norm(v); cls=(n@v)%3; ul=u@e
        m=[ul[cls==c].mean() for c in range(3)]; amp=np.ptp(m)
        if best is None or amp>best[0]: best=(amp,v,cls,ul)
    return best
a=ax[1,0]
for j,(tag,c,lab) in enumerate([('Zr50Ru2Hf2_c1','#2a78d6','Zr$_{50}$Ru$_2$Hf$_2$'),('Zr50Ru4_c1','#0b0b0b','base (config. 1)')]):
    amp,v,cls,ul=plane(tag)
    # order classes so that the 'fixed' plane (smallest |mean|) is first
    means=[ul[cls==k].mean() for k in range(3)]; order=np.argsort(np.abs(means))
    for i,k in enumerate(order):
        y=ul[cls==k]; xx=i+(j-0.5)*0.32+np.random.default_rng(1).uniform(-0.07,0.07,len(y))
        a.plot(xx,y,'o',ms=3,mfc=c,mec='white',mew=0.3,alpha=0.85,label=lab if i==0 else None)
d0=np.sqrt(3)/12*3.55
for yv in (d0,-d0): a.axhline(yv,color='#9a9890',ls='--',lw=0.8)
a.text(2.45,d0+0.02,'ideal ω',fontsize=6,color='#7a7974',ha='right')
a.set_xticks([0,1,2]); a.set_xticklabels(['plane A\n(fixed in ω)','plane B','plane C']); a.set_ylabel('displacement along ⟨111⟩ (Å)'); a.legend(fontsize=6.5,loc='lower left')
a.axhline(0,color='#c3c2b7',lw=0.6)
a=ax[1,1]
V=U.sort_values('omega_frac_max')
a.barh(range(len(V)),V.omega_frac_max*100,zorder=2,color=[FAM[f][0] for f in V.fam_true],height=0.75,edgecolor='white',lw=0.5)
a.axvspan(0,7.5,color='#e8e6df',lw=0,zorder=0); a.axvline(2.63,color='#0b0b0b',ls='--',lw=0.9,zorder=1); a.text(42,4,'random displacement field:\nmean 2.6% (dashed),\n99th percentile 7.5% (shaded)',fontsize=6)
a.axvline(40,color='#9a9890',ls=':',lw=0.9)
a.set_yticks(range(len(V))); a.set_yticklabels([s+(' †' if st=='B off-site' else '') for s,st in zip(V.short,V.state)],fontsize=6)
a.set_xlabel('share of static displacement in one ω variant (%)')
for i,(axx,l) in enumerate(zip(ax.ravel(),'abcd')): axx.text(-0.14,1.03,f'({l})',transform=axx.transAxes,fontweight='bold')
fig.tight_layout(); fig.savefig(OUT+'Fig2_relaxation_omega.png'); fig.savefig(OUT+'Fig2_relaxation_omega.pdf'); plt.close()

# ---------- Fig 3: C' and G/B vs VEC ----------
nb=U[U.state=='near-bcc']
fig,ax=plt.subplots(1,2,figsize=(7.1,3.0))
for a,y,yl in [(ax[0],'Cp',"$C'=(C_{11}-C_{12})/2$ (GPa)"),(ax[1],'G_over_B','Pugh ratio $G/B$')]:
    s=stats.linregress(nb.VEC,nb[y]); xx=np.linspace(4.05,4.43,50)
    n=len(nb); xm=nb.VEC.mean(); sxx=((nb.VEC-xm)**2).sum(); res=nb[y]-(s.intercept+s.slope*nb.VEC); se=np.sqrt((res**2).sum()/(n-2))
    band=stats.t.ppf(0.975,n-2)*se*np.sqrt(1/n+(xx-xm)**2/sxx)
    a.fill_between(xx,s.intercept+s.slope*xx-band,s.intercept+s.slope*xx+band,color='#e8e6df',lw=0)
    a.plot(xx,s.intercept+s.slope*xx,color='#7a7974',lw=1)
    scatter(a,'VEC',y,offs={'4Y':(3,-7),'4Sc':(-14,2),'Hf':(3,2)})
    a.set_xlabel('valence-electron concentration (e/atom)'); a.set_ylabel(yl)
    a.text(0.97 if y!='Cp' else 0.03,0.95,f'near-BCC fit: r = {s.rvalue:.2f}, slope = {s.slope:.{0 if y=="Cp" else 2}f}',transform=a.transAxes,fontsize=6.5,va='top',ha='right' if y!='Cp' else 'left')
famlegend(ax[0],loc='lower right')
for axx,l in zip(ax,'ab'): axx.text(-0.15,1.03,f'({l})',transform=axx.transAxes,fontweight='bold')
fig.tight_layout(); fig.savefig(OUT+'Fig3_Cprime_GB_VEC.png'); fig.savefig(OUT+'Fig3_Cprime_GB_VEC.pdf'); plt.close()

# ---------- Fig 4: electronic ----------
fig,ax=plt.subplots(2,2,figsize=(7.1,5.4))
W=U.dropna(subset=['NEF_tot'])
scatter(ax[0,0],'VEC','NEF_tot',data=W,offs={'4Y':(3,-7)}); ax[0,0].set_xlabel('valence-electron concentration (e/atom)'); ax[0,0].set_ylabel('$N(E_F)$ (states eV$^{-1}$ atom$^{-1}$)')
scatter(ax[0,1],'NEF_tot','Cp',data=W); ax[0,1].set_xlabel('$N(E_F)$ (states eV$^{-1}$ atom$^{-1}$)'); ax[0,1].set_ylabel("$C'$ (GPa)")
s=stats.linregress(W[W.state=='near-bcc'].NEF_tot,W[W.state=='near-bcc'].Cp); ax[0,1].text(0.97,0.95,f'near-BCC: r = {s.rvalue:.2f}',transform=ax[0,1].transAxes,ha='right',va='top',fontsize=6.5)
scatter(ax[1,0],'NEF_tot','G_over_B',data=W); ax[1,0].set_xlabel('$N(E_F)$ (states eV$^{-1}$ atom$^{-1}$)'); ax[1,0].set_ylabel('$G/B$')
s=stats.linregress(W.NEF_tot,W.G_over_B); ax[1,0].text(0.97,0.95,f'all: r = {s.rvalue:.2f}',transform=ax[1,0].transAxes,ha='right',va='top',fontsize=6.5)
famlegend(ax[1,0],loc='lower left')
a=ax[1,1]
def sm(y,dE,s=0.10):
    k=np.arange(-4*s,4*s+dE,dE); g=np.exp(-k**2/2/s**2); g/=g.sum(); return np.convolve(y,g,'same')
for tag,c,lab in [('Zr50Ru4_c1','#0b0b0b','base (config. 1)'),('Zr48Ru4Re2_c1','#eb6834','Zr$_{48}$Ru$_4$Re$_2$ (near-BCC)'),
                   ('Zr50Ru2Hf2_c1','#2a78d6','Zr$_{50}$Ru$_2$Hf$_2$ (ω-collapsed)'),('Zr48Ru2Y4_c1','#1baf7a','Zr$_{48}$Ru$_2$Y$_4$ (distorted)')]:
    d=np.load(f'data/dos/{tag}.npz'); E=d['E']; dE=E[1]-E[0]
    a.plot(E,sm(d['tot'],dE)/54,color=c,lw=1.1,label=lab)
a.axvline(0,color='#9a9890',lw=0.8,ls='--'); a.set_xlim(-4.5,0.9); a.set_ylim(0,2.1); a.set_xlabel('$E-E_F$ (eV)'); a.set_ylabel('total DOS (states eV$^{-1}$ atom$^{-1}$)')
a.legend(fontsize=6.3,loc='upper left'); 
for axx,l in zip(ax.ravel(),'abcd'): axx.text(-0.15,1.03,f'({l})',transform=axx.transAxes,fontweight='bold')
fig.tight_layout(); fig.savefig(OUT+'Fig5_electronic.png'); fig.savefig(OUT+'Fig5_electronic.pdf'); plt.close()

# ---------- Fig 5: local environment, Pugh-Pettifor & distortion ----------
fig,ax=plt.subplots(1,3,figsize=(7.1,3.0),gridspec_kw={'width_ratios':[0.8,1,1]})
E=pd.read_csv('results/zr_local_environment.csv')
a=ax[0]; g=E.groupby('nRu').u
for k,(nr,vals) in enumerate(g):
    q=np.percentile(vals,[25,50,75]); a.plot([nr,nr],[q[0],q[2]],color='#2a78d6',lw=3,alpha=0.5,solid_capstyle='butt')
    a.plot(nr,vals.mean(),'o',color='#2a78d6',ms=5); a.text(nr,0.56,f'n={len(vals)}',ha='center',fontsize=5.8,color='#52514e')
a.set_xticks([0,1,2,3]); a.set_xlabel('Ru atoms in 1st+2nd shell of Zr'); a.set_ylabel('Zr static displacement (Å)'); a.set_ylim(0,0.6)
a=ax[1]; scatter(a,'G_over_B','Cauchy',fs=5.8,offs={'4Ta':(4,-3),'Mo':(3,4),'Nb':(5,-1),'4Y':(-12,-8),'Ti':(-9,4)})
a.set_xlabel('$G/B$'); a.set_ylabel('$C_{12}-C_{44}$ (GPa)')
a.set_xlim(0.295,0.42); a.text(0.415,49.5,'all: $G/B$ < 0.57,\n$C_{12}-C_{44}$ > 0',ha='right',va='top',fontsize=6,color='#7a7974')
a=ax[2]; scatter(a,'u_rms','G_over_B',fs=5.8,offs={'4Y':(3,-7)}); a.set_xlabel('RMS static displacement $\\bar u$ (Å)'); a.set_ylabel('$G/B$')
from matplotlib.lines import Line2D
h=[Line2D([],[],marker=m,ls='none',mfc=c,mec=c,ms=s*0.8,label=l) for f,(c,m,s,l) in FAM.items()]+[Line2D([],[],marker='o',ls='none',mfc='white',mec='#52514e',mew=1.3,ms=5,label='ω-collapsed / strongly distorted')]
fig.legend(handles=h,loc='lower center',ncol=5,fontsize=6.3,frameon=False)
for axx,l in zip(ax,'abc'): axx.text(-0.2,1.03,f'({l})',transform=axx.transAxes,fontweight='bold')
fig.tight_layout(rect=(0,0.07,1,1)); fig.savefig(OUT+'Fig4_local_pugh_distortion.png'); fig.savefig(OUT+'Fig4_local_pugh_distortion.pdf'); plt.close()

# ---------- Fig 6: design map ----------
fig,ax=plt.subplots(figsize=(4.6,3.6))
base=U[U.fam_true=='BASE']; bC=base.Cp.mean(); bu=base.u_rms.mean()
ax.fill_between([0.25,base.u_rms.min()],base.Cp.max(),40,color='#f6d9c6',lw=0)
ax.text(0.2515,36.0,'stiffer $C\'$ and smaller static\ndistortion than both\nZr$_{50}$Ru$_4$ configurations',fontsize=6.3,color='#9b3d12')
from matplotlib.patches import Rectangle
ax.add_patch(Rectangle((base.u_rms.min(),base.Cp.min()),np.ptp(base.u_rms),np.ptp(base.Cp),fill=False,ec='#7a7974',lw=0.8,ls='--'))
scatter(ax,'u_rms','Cp',fs=6.5,offs={'Re':(3,3),'4Y':(4,-8),'Ti':(3,-7),'Sc':(-9,4),'V':(4,-7)},skip={('BASE','base'),('Z','Hf'),('Z','Al')})
ax.annotate('base (2 configurations)\nand Hf (Z)',(0.312,26.6),xytext=(0.335,26.2),fontsize=6,color='#3a3a38',arrowprops=dict(arrowstyle='-',color='#9a9890',lw=0.6))
ax.annotate('Al (Z)',(0.321,29.2),xytext=(0.333,30.5),fontsize=6,color='#3a3a38',arrowprops=dict(arrowstyle='-',color='#9a9890',lw=0.6))
ax.set_xlabel('RMS static displacement $\\bar u$ (Å)'); ax.set_ylabel("$C'$ (GPa)"); ax.set_xlim(0.25,0.50); ax.set_ylim(18,38)
famlegend(ax,loc='lower right')
fig.tight_layout(); fig.savefig(OUT+'Fig6_design_map.png'); fig.savefig(OUT+'Fig6_design_map.pdf'); plt.close()

# ---------------- Supplementary figures
plt.rcParams.update({'font.size':7})
# S1: PDOS grid
reps=U[U.id.map(lambda i: DS[i]['dos_file'] is not None)].sort_values(['fam_true','VEC'])
n=len(reps); nc=4; nr=int(np.ceil(n/nc))
fig,ax=plt.subplots(nr,nc,figsize=(7.1,1.45*nr),sharex=True)
for a,(_,r) in zip(ax.ravel(),reps.iterrows()):
    e=DS[r.id]; d=dict(np.load(e['dos_file'])); sp=e['species']; cnt=dict(zip(sp,e['counts']))
    E=d['E']; dE=E[1]-E[0]
    a.fill_between(E,sm(d['tot'],dE)/sum(cnt.values()),color='#e8e6df',lw=0,label='total/atom')
    a.plot(E,sm(d['Zr_t2g']+d['Zr_eg'],dE)/cnt['Zr'],color='#0b0b0b',lw=0.8,label='Zr-d/atom')
    a.plot(E,sm(d['Ru_t2g']+d['Ru_eg'],dE)/cnt['Ru'],color='#2a78d6',lw=0.8,label='Ru-d/atom')
    X=[e for e in sp if e not in ('Zr','Ru')]
    if X:
        x=X[0]; a.plot(E,sm(d[x+'_t2g']+d[x+'_eg']+d[x+'_s']+d[x+'_p'],dE)/cnt[x],color='#eb6834',lw=0.8,label='X (s+p+d)/atom')
    a.axvline(0,color='#9a9890',lw=0.6,ls='--'); a.set_xlim(-6,0.9); a.set_ylim(0,2.6)
    a.set_title(f"{r.label} ({r.state.replace('omega','ω')})",fontsize=6.5)
for a in ax.ravel()[n:]: a.axis('off')
ax[0,0].legend(fontsize=5.5,loc='upper left')
for a in ax[-1]: a.set_xlabel('$E-E_F$ (eV)')
fig.tight_layout(); fig.savefig(OUT+'FigS1_PDOS_all.png'); plt.close()
# S2: trajectories
T=pd.read_csv('results/trajectories.csv').rename(columns={'id':'tag','energy_eV':'E'})
fig,ax=plt.subplots(1,2,figsize=(7.1,2.8))
for tag,g in T.groupby('tag'):
    if g.step.max()<10: continue
    g=g.sort_values('step'); ax[0].plot(g.step,(g.E-g.E.iloc[-1])*1000/54,lw=0.6,color='#7a7974',alpha=0.7)
    ax[1].plot(g.u_rms,g.omega_amp,lw=0.6,color='#7a7974',alpha=0.7); ax[1].plot(g.u_rms.iloc[-1],g.omega_amp.iloc[-1],'o',ms=2.5,color='#2a78d6')
ax[0].set_yscale('symlog',linthresh=1); ax[0].set_xlabel('ionic step'); ax[0].set_ylabel('$E-E_\\mathrm{final}$ (meV/atom)')
ax[1].set_xlabel('$\\bar u$ (Å)'); ax[1].set_ylabel('$u_\\omega$ (Å)')
fig.tight_layout(); fig.savefig(OUT+'FigS2_trajectories_all.png'); plt.close()
print('figures written to figures/')

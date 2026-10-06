"""Publication figures from saved coordinates only; no embedding/model fit."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[2]; D=P/'runs/zymo_fecal_20261005/derived'; O=P/'docs/research/figures/zymo_20261005'; O.mkdir(parents=True,exist_ok=True)
f=pd.read_csv(D/'joined_read_metadata.tsv.gz',sep='\t'); order=np.random.default_rng(42).permutation(len(f))
plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.labelsize':9,'figure.facecolor':'white'})
species=f.KC_species.value_counts().drop('unresolved').head(8).index.tolist(); colors=list(plt.get_cmap('tab10').colors)
for kind in ['taxa','length','gc']:
 fig,axes=plt.subplots(2,2,figsize=(12,8.8),layout='constrained')
 for row,method in enumerate(['tsne','umap']):
  for col,v in enumerate([1,2]):
   ax=axes[row,col]; x=f[f'{method}{v}_x'].values; y=f[f'{method}{v}_y'].values
   if kind=='taxa':
    c=np.array([colors[species.index(t)] if t in species else (.76,.76,.76) for t in f.KC_species]); ax.scatter(x[order],y[order],c=c[order],s=.35,alpha=.6,linewidths=0,rasterized=True)
    for j,t in enumerate(species):
     g=f[f.KC_species==t]; cx=g[f'{method}{v}_x'].median(); cy=g[f'{method}{v}_y'].median(); offset=(-32,10) if j==6 else ((24,-17) if j==7 else (0,0))
     ax.annotate(str(j+1),(cx,cy),xytext=offset,textcoords='offset points',ha='center',va='center',fontsize=10,weight='bold',bbox={'facecolor':'white','alpha':.9,'edgecolor':'none','pad':1.5},arrowprops={'arrowstyle':'-','color':'#555555','lw':.7} if j>=6 else None)
   else:
    c=np.log10(f.read_len.to_numpy()) if kind=='length' else f.gc_exact.to_numpy(); lo,hi=(3,5) if kind=='length' else (.15,.75)
    sc=ax.scatter(x[order],y[order],c=c[order],cmap='viridis' if kind=='length' else 'coolwarm',vmin=lo,vmax=hi,s=.35,alpha=.7,linewidths=0,rasterized=True)
   ax.set_title(f'{"t-SNE" if method=="tsne" else "UMAP"} | version {v}');ax.set_xticks([]);ax.set_yticks([])
   ax.set_xlabel('Display dimension 1');ax.set_ylabel('Display dimension 2')
 if kind!='taxa': fig.colorbar(sc,ax=axes.ravel().tolist(),shrink=.78,label='log10 length (bp); 3 = 1 kb, 4 = 10 kb, 5 = 100 kb' if kind=='length' else 'GC fraction (fixed color scale)')
 else:
  from matplotlib.lines import Line2D
  handles=[Line2D([],[],marker='o',linestyle='',color=colors[i],label=f'{i+1}. {t}') for i,t in enumerate(species)]
  fig.legend(handles=handles,loc='outside lower center',ncol=2,fontsize=8,frameon=False)
 fig.savefig(O/f'{kind}_comparison.png',dpi=160); plt.close(fig)
# Preserve the prior exact bin overlay as the report's bin evidence figure.
import shutil
shutil.copy2(D/'comparison_bins.png',O/'bin_comparison.png')
# Compact supplemental observed stats; no read IDs exported.
extra={'length_strata':[], 'rare_target_calls':{},'candidate_bin_labels':__import__('json').loads((D/'targeted_bin_labels.json').read_text())}
for lo,hi in [(0,2000),(2000,5000),(5000,10000),(10000,30000),(30000,100000)]:
 q=f[(f.read_len>lo)&(f.read_len<=hi)];extra['length_strata'].append({'lower_exclusive_bp':lo,'upper_inclusive_bp':hi,'reads':len(q),'KC_species_resolved_fraction':float(q.KC_species.ne('unresolved').mean())})
for target in ['Porphyromonas gingivalis','Peptostreptococcus','Gemella']:
 extra['rare_target_calls'][target]={tool:{'total':int(f[tool+'_species'].str.contains(target,regex=False).sum()),'labels':f.loc[f[tool+'_species'].str.contains(target,regex=False),tool+'_species'].value_counts().to_dict()} for tool in ['M','K','C']}
E=P/'docs/computing/evidence/zymo_20261005'; (E/'supplemental_findings.json').write_text(__import__('json').dumps(extra,indent=2),encoding='utf8')
print('Four figures and supplemental aggregate evidence prepared.')
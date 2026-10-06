"""Join by validated read identity, audit raw bins, and render diagnostic views."""
from pathlib import Path
import json,pickle,io
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
P=Path(__file__).resolve().parents[2]/'runs/zymo_fecal_20261005'; D=P/'derived'
f=pd.read_csv(D/'plot_read_metadata.tsv.gz',sep='\t',index_col='read_id')
m=pd.read_csv(D/'read_to_lrb.tsv',sep='\t',index_col='read_id')
assert m.index.is_unique and set(f.index)==set(m.index)
m=m.loc[f.index]; assert np.array_equal(f.read_len,m.read_len)
assert np.max(abs(f.gc-m.gc_exact))<=0.00500001
f=f.join(m.drop(columns='read_len')); f.to_csv(D/'joined_read_metadata.tsv.gz',sep='\t')
class BasicUnpickler(pickle.Unpickler):
 def find_class(self,*args): raise ValueError('Pickle globals forbidden')
p=BasicUnpickler(io.BytesIO((P/'raw/zymo_lrb/binning_result.pkl').read_bytes())).load()
bins=np.loadtxt(P/'raw/zymo_lrb/bins.txt',dtype=np.int16)
lens=np.loadtxt(P/'raw/zymo_lrb/lengths.txt',dtype=np.int64)
a=np.full(len(bins),-1,np.int16); duplicates=0
for k,v in p.items():
 ix=np.asarray(v); assert len(np.unique(ix))==len(ix); duplicates+=int((a[ix]>=0).sum()); a[ix]=k
assert not duplicates and (a>=0).all() and np.array_equal(a,bins)
assert np.array_equal(bins[f.full_index.values],f.lrb_bin.values)
summary={'full_reads':len(bins),'full_bases':int(lens.sum()),'plot_read_fraction':len(f)/len(bins),'plot_base_fraction':f.read_len.sum()/lens.sum(),'pickle_bins_match_bins_txt':True,'duplicate_bin_memberships':duplicates,'unbinned':int((a<0).sum()),'first_plotted_full_index':int(f.full_index.min()),'last_plotted_full_index':int(f.full_index.max()),'source_read_fraction_before_last_selected':(f.full_index.max()+1)/len(bins)}
records=[]
for b,group in f.groupby('lrb_bin'):
 vc=group.KC_species.value_counts(); known=group[group.KC_species!='unresolved']; top=known.KC_species.value_counts(); ng=len(group)
 records.append({'bin':b,'full_reads':int((bins==b).sum()),'full_read_pct':100*np.mean(bins==b),'plot_reads':ng,'plot_read_pct':100*ng/len(f),'sampling_pct':100*ng/(bins==b).sum(),'median_len':group.read_len.median(),'mean_gc':group.gc_exact.mean(),'KC_species_resolved_pct':100*len(known)/ng,'top_KC_species':top.index[0] if len(top) else 'none','top_species_pct_all_bin_reads':100*top.iloc[0]/ng if len(top) else 0,'top3_KC_species':'; '.join(f'{t}: {n} ({100*n/ng:.1f}%)' for t,n in top.head(3).items())})
bt=pd.DataFrame(records).sort_values('full_reads',ascending=False); bt.to_csv(D/'bin_summary.csv',index=False)
sampling=[]
for i in range(10):
 lo=i*len(bins)//10; hi=(i+1)*len(bins)//10; selected=f[(f.full_index>=lo)&(f.full_index<hi)]
 sampling.append({'decile':i+1,'full_start_zero_based':lo,'full_end_exclusive':hi,'full_reads':hi-lo,'selected_reads':len(selected),'selection_pct':100*len(selected)/(hi-lo),'full_median_len':float(np.median(lens[lo:hi]))})
pd.DataFrame(sampling).to_csv(D/'sampling_by_source_decile.csv',index=False); summary['sampling_by_decile']=sampling
idx=np.sort(f.full_index.to_numpy()); dif=np.diff(idx); summary['contiguous_selected_runs']=int(1+(dif>1).sum()); summary['adjacent_selected_pairs']=int((dif==1).sum())
rng=np.random.default_rng(42); qidx=rng.choice(len(f),10000,replace=False); nbr={}; metrics={}
for v in [1,2]:
 for method in ['tsne','umap']:
  name=f'{method}{v}'; xy=f[[name+'_x',name+'_y']].to_numpy(); _,ind=cKDTree(xy).query(xy[qidx],k=17,workers=2); ind=np.array([row[row!=qi][:15] for qi,row in zip(qidx,ind)]); nbr[name]=ind
  same=(f.lrb_bin.values[qidx,None]==f.lrb_bin.values[ind]); lab=f.KC_species.to_numpy(); valid=(lab[qidx,None]!='unresolved')&(lab[ind]!='unresolved')
  metrics[name]={'same_lrb_bin_neighbor_fraction':float(same.mean()),'same_KC_species_fraction_among_both_resolved':float((lab[qidx,None]==lab[ind])[valid].mean()),'resolved_neighbor_pairs':int(valid.sum()),'median_abs_log2_length_difference':float(np.median(abs(np.log2(f.read_len.values[qidx,None]/f.read_len.values[ind])))),'median_abs_gc_difference':float(np.median(abs(f.gc_exact.values[qidx,None]-f.gc_exact.values[ind])))}
for a,b in [('tsne1','tsne2'),('umap1','umap2'),('tsne1','umap1'),('tsne2','umap2')]: metrics[a+'_'+b+'_neighbor_overlap']=float(np.mean([len(set(x)&set(y))/15 for x,y in zip(nbr[a],nbr[b])]))
summary['neighbor_diagnostics_k15_10000_queries_seed42']=metrics
(D/'comparison_summary.json').write_text(json.dumps(summary,indent=2))
colors=list(plt.get_cmap('tab20').colors)+list(plt.get_cmap('Set2').colors)
# All points retained, draw order deterministically shuffled to reduce categorical occlusion.
order=rng.permutation(len(f)); species=f.KC_species.value_counts().drop('unresolved').head(8).index.tolist()
for kind in ['bins','taxa','gradients']:
 fig,axes=plt.subplots(2,2,figsize=(13,10),layout='constrained')
 for row,method in enumerate(['tsne','umap']):
  for col,v in enumerate([1,2]):
   ax=axes[row,col]; name=f'{method}{v}'; x=f[name+'_x'].values; y=f[name+'_y'].values
   if kind=='bins':
    c=np.array([colors[int(b)%len(colors)] for b in f.lrb_bin]); ax.scatter(x[order],y[order],c=c[order],s=.3,alpha=.5,rasterized=True,linewidths=0)
    for b,g in f.groupby('lrb_bin'):
     if len(g)>=150: ax.text(g[name+'_x'].median(),g[name+'_y'].median(),str(b),fontsize=8,weight='bold',bbox={'facecolor':'white','alpha':.7,'edgecolor':'none','pad':1})
   elif kind=='taxa':
    ids=np.array([species.index(t) if t in species else -1 for t in f.KC_species]); c=np.array([colors[i] if i>=0 else (.75,.75,.75) for i in ids]); ax.scatter(x[order],y[order],c=c[order],s=.3,alpha=.5,rasterized=True,linewidths=0)
    for i,t in enumerate(species):
     g=f[f.KC_species==t]; ax.text(g[name+'_x'].median(),g[name+'_y'].median(),str(i+1),fontsize=9,weight='bold',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
   else:
    c=np.log10(f.read_len.values) if col==0 else f.gc_exact.values; sc=ax.scatter(x[order],y[order],c=c[order],cmap='viridis' if col==0 else 'coolwarm',s=.3,alpha=.6,rasterized=True,linewidths=0); fig.colorbar(sc,ax=ax,label='log10 read length (bp)' if col==0 else 'GC fraction')
   ax.set_title(f'{method.upper()} version {v}'); ax.set_xlabel('Display dimension 1'); ax.set_ylabel('Display dimension 2'); ax.set_xticks([]); ax.set_yticks([])
 title={'bins':'LRBinner full-run bins on TaxoViz coordinates','taxa':'KrakenUniq–Centrifuger species agreement (diagnostic labels)','gradients':'Length / GC diagnostic gradients'}[kind]
 fig.suptitle(title+' | 290,227 identical reads',fontsize=14)
 if kind=='taxa': fig.supxlabel('\n'.join([',  '.join(f'{i+1} {t}' for i,t in enumerate(species[:4])),',  '.join(f'{i+5} {t}' for i,t in enumerate(species[4:]))])+'\nGray: other taxa or unresolved',fontsize=9)
 fig.savefig(D/f'comparison_{kind}.png',dpi=160); plt.close(fig)
print(bt.to_string(index=False)); print(json.dumps(summary,indent=2))
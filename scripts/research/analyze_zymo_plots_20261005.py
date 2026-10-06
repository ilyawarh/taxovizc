"""Reproduce descriptive review from saved Plotly traces. No embedding/model is fitted."""
from pathlib import Path
import json,base64,hashlib,pickle,collections,re,sys
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parents[2]/'runs/zymo_fecal_20261005'
D=P/'derived'; D.mkdir(exist_ok=True)
def array(v):
 return np.frombuffer(base64.b64decode(v['bdata']),dtype=v['dtype']) if isinstance(v,dict) else np.asarray(v)
def parse(v):
 s=(P/f'raw/zymofec_{v}.html').read_text(encoding='utf8'); i=s.index('[',s.rfind('Plotly.newPlot(')); traces,_=json.JSONDecoder().raw_decode(s[i:]); del s
 frames={}
 for axis in ['x','x2','x3','x4']:
  parts=[]
  for t in traces:
   if t.get('xaxis')!=axis: continue
   rows=[dict(z.split(': ',1) for z in txt.split('<br>')) for txt in t['text']]
   q=pd.DataFrame(rows).rename(columns={'read':'read_id','len':'read_len'}); q['x']=array(t['x']); q['y']=array(t['y']); q['category']=t['legendgroup']; parts.append(q)
  frames[axis]=pd.concat(parts,ignore_index=True).set_index('read_id').sort_index()
  assert frames[axis].index.is_unique
 assert all(frames[a].index.equals(frames['x'].index) for a in frames)
 for a,b in [('x','x3'),('x2','x4')]: assert frames[a][['read_len','gc','M','K','C','category']].equals(frames[b][['read_len','gc','M','K','C','category']])
 for a,b in [('x','x2'),('x3','x4')]: assert np.array_equal(frames[a][['x','y']].values,frames[b][['x','y']].values)
 f=frames['x'][['read_len','gc']].astype({'read_len':int,'gc':float})
 for rank,ax in [('genus','x'),('species','x2')]:
  for c in ['M','K','C','category']: f[c+'_'+rank]=frames[ax][c]
 for method,ax in [('tsne','x'),('umap','x3')]:
  for c in ['x','y']: f[f'{method}{v}_{c}']=frames[ax][c]
 return f
f=parse(1); g=parse(2)
(D/'wanted_read_ids.txt').write_text('\n'.join(f.index)+'\n')
cols=[c for c in f if not c.startswith(('tsne','umap'))]
assert f.index.equals(g.index) and f[cols].equals(g[cols])
f=f.join(g[[c for c in g if c.startswith(('tsne','umap'))]])
summary={'n_reads':len(f),'duplicate_plot_ids':0,'both_versions_identical_metadata':True,'total_bases':int(f.read_len.sum()),'read_length_quantiles':f.read_len.quantile([0,.1,.25,.5,.75,.9,.99,1]).to_dict(),'mean_gc_rounded':f.gc.mean(),'agreement':{},'actual_availability':{},'hidden_conflicts':{}}
for rank in ['genus','species']:
 q=f[[x+'_'+rank for x in ['M','K','C']]]; present=q.ne('unclassified'); n=present.sum(axis=1)
 summary['agreement'][rank]=f['category_'+rank].value_counts().to_dict()
 summary['actual_availability'][rank]={x:int(present[x+'_'+rank].sum()) for x in ['M','K','C']}
 conflicts=q.apply(lambda x:len(set(x)-{'unclassified'})>1,axis=1)
 summary['hidden_conflicts'][rank]={c:{'total':int((f['category_'+rank]==c).sum()),'with_conflict':int(((f['category_'+rank]==c)&conflicts).sum())} for c in ['M_only','K_only','C_only','KC','MK','MC','MKC','none']}
 f['KC_'+rank]=np.where((f['K_'+rank]==f['C_'+rank]) & f['K_'+rank].ne('unclassified'), f['K_'+rank], 'unresolved')
 for tool in ['M','K','C','KC']:
  c=tool+'_'+rank
  t=f.groupby(c).agg(reads=('read_len','size'),bases=('read_len','sum'),median_len=('read_len','median'),mean_gc=('gc','mean')).sort_values('reads',ascending=False)
  t['read_pct']=100*t.reads/len(f); t['base_pct']=100*t.bases/f.read_len.sum(); t.to_csv(D/f'taxa_{c}.csv')
f.to_csv(D/'plot_read_metadata.tsv.gz',sep='\t')
(D/'plot_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
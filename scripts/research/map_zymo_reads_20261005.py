"""Read-only SSH extraction of plotted IDs in original FASTQ order; no remote files written."""
import pathlib, subprocess, shlex
p=pathlib.Path('/mnt/c/Users/1/Desktop/SMRL/runs/zymo_fecal_20261005')
remote=r'''
import sys, json, time, os
from itertools import zip_longest
base='/mnt/raid0/Toxins/test/lrbinner/zymo_fecal_run/'
wanted=set(sys.stdin.buffer.read().splitlines())
found=set(); dup=0; n=0; bases=0; matched_bases=0; start=time.time()
fq=base+'Zymo_fecal.trimmed.fastq'
print('read_id\tfull_index\tlrb_bin\tread_len\tgc_exact')
with open(fq,'rb',buffering=8*1024*1024) as f, open(base+'zymo_lrb/bins.txt') as bins, open(base+'zymo_lrb/lengths.txt') as lens:
 while True:
  h=f.readline()
  if not h: break
  seq=f.readline().rstrip(b'\r\n'); plus=f.readline(); qual=f.readline().rstrip(b'\r\n')
  if not h.startswith(b'@') or not plus.startswith(b'+') or len(seq)!=len(qual): raise ValueError(('Malformed four-line FASTQ',n))
  b=bins.readline().strip(); stored=lens.readline().strip()
  if not b or int(stored)!=len(seq): raise ValueError(('LRB length/cardinality mismatch',n))
  rid=h[1:].split()[0]; bases+=len(seq)
  if rid in wanted:
   dup+=rid in found; found.add(rid); matched_bases+=len(seq)
   seq=seq.upper(); gc=(seq.count(b'G')+seq.count(b'C'))/len(seq)
   print(rid.decode(),n,b,len(seq),gc,sep='\t')
  n+=1
  if n%1000000==0: print(json.dumps({'processed':n,'matched':len(found),'elapsed_s':round(time.time()-start)}),file=sys.stderr,flush=True)
 if bins.readline() or lens.readline(): raise ValueError('Extra LRB rows')
summary={'full_reads':n,'full_bases':bases,'wanted':len(wanted),'matched':len(found),'duplicate_matched_ids':dup,'missing':len(wanted-found),'matched_bases':matched_bases,'elapsed_s':time.time()-start,'full_fastq_bytes':os.path.getsize(fq),'full_fastq_mtime':os.path.getmtime(fq),'validation':'all records strict four-line and all lengths match LRB lengths.txt; IDs unique among matched reads; global nonmatched duplicate IDs not tested'}
print(json.dumps(summary),file=sys.stderr,flush=True)
'''
with (p/'derived/read_to_lrb.tsv').open('wb') as out, (p/'derived/read_mapping_log.jsonl').open('wb') as err:
 r=subprocess.run(['ssh','-o','BatchMode=yes','dgx03-ngate','python3 -c '+shlex.quote(remote)],input=(p/'derived/wanted_read_ids.txt').read_bytes(),stdout=out,stderr=err)
raise SystemExit(r.returncode)
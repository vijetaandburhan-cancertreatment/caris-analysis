#!/usr/bin/env python3
"""Independent plain-Python FASTQ/str.find recount on extracted candidate pairs."""
import argparse,collections,gzip,hashlib,json,pathlib,time

def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def nm(s):return s[:-2] if s.endswith(('/1','/2')) else s

def fastq(path):
 with gzip.open(path,'rt') as f:
  while True:
   h=f.readline()
   if not h:return
   s=f.readline().strip();plus=f.readline();q=f.readline().strip();assert h.startswith('@') and plus.startswith('+') and len(s)==len(q)
   yield nm(h[1:].split()[0]),s,q

def main():
 p=argparse.ArgumentParser();p.add_argument('--audit',required=True);a=p.parse_args();root=pathlib.Path(a.audit);primary=json.load(open(root/'junction-read-audit.json'));rows=primary['rows'];t0=time.time();words=collections.defaultdict(list);named={r['audit_key']:{nm(n) for n in r['read_names']} for r in rows}
 for r in rows:
  for ptn in r['junction_patterns']:
   for w in {ptn['sequence'],rc(ptn['sequence'])}:words[w].append((r['audit_key'],ptn['arm_bases']))
 counts=collections.Counter();families=collections.defaultdict(set);n=0
 for r1,r2 in zip(fastq(root/'audit-selected.R1.fastq.gz'),fastq(root/'audit-selected.R2.fastq.gz'),strict=True):
  assert r1[0]==r2[0];n+=1;good=set()
  for rn,s,q in [r1,r2]:
   for w,targets in words.items():
    start=s.find(w)
    while start>=0:
     quality=min(ord(c)-33 for c in q[start:start+len(w)])
     for t,arm in targets:
      for threshold in [20,30]:
       if quality>=threshold:good.add((t,str(arm),f'Q{threshold}'))
     start=s.find(w,start+1)
  sh=hashlib.sha256((r1[1]+'|'+r2[1]).encode()).hexdigest()
  for key in good:
   counts[key+('pairs',)]+=1;counts[key+('caller_named_pairs',)]+=r1[0] in named[key[0]];families[key].add(sh)
 checks=[]
 for r in rows:
  t=r['audit_key']
  for arm in ['15','20','25']:
   for q in ['Q20','Q30']:
    s=primary['statistics'][t]['raw_exact'][arm][q];key=(t,arm,q)
    for field in ['pairs','caller_named_pairs','sequence_families']:
     actual=len(families[key]) if field=='sequence_families' else counts[key+(field,)];checks.append({'candidate':t,'arm':int(arm),'quality':q,'field':field,'expected':s[field],'recounted':actual,'match':actual==s[field]})
 out={'status':'passed' if all(x['match'] for x in checks) else 'failed','candidate_pairs_recounted':n,'checks':checks,'elapsed_seconds':time.time()-t0,'method':'Separate gzip four-line FASTQ parser plus exhaustive str.find of each exact forward/reverse marker. Reproduces Q20/Q30 pair counts, caller-name subset and full read-pair sequence families. No pysam/Aho-Corasick used.','limitations':['Uses the exhaustive primary full-scan candidate extraction, so this recount alone is not a second search for omitted-negative raw pairs. Public matcher regression separately tests the full-scan matching implementation.']};(root/'independent-recount.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:out[k] for k in ['status','candidate_pairs_recounted','elapsed_seconds']}));assert out['status']=='passed'
if __name__=='__main__':main()

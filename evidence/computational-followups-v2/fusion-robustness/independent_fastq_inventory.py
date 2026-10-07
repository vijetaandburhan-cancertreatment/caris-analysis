#!/usr/bin/env python3
"""Modest-set FASTQ QA, separate stdlib parser; no inference of original molecules."""
import argparse,collections,gzip,hashlib,itertools,json,pathlib,time

def records(path):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt') as f:
  while True:
   h=f.readline()
   if not h:return
   s=f.readline().rstrip('\r\n');plus=f.readline();q=f.readline().rstrip('\r\n')
   if not h.startswith('@') or not plus.startswith('+') or len(s)!=len(q) or not q:raise ValueError('FASTQ header/plus/sequence-quality failure')
   n=h[1:].split()[0];n=n[:-2] if n.endswith(('/1','/2')) else n
   yield n,s,q

def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def checksum(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--r1',required=True);p.add_argument('--r2',required=True);p.add_argument('--out',required=True);a=p.parse_args();t=time.time();names=collections.Counter();families=collections.Counter();normalized=collections.Counter();lengths=[collections.Counter(),collections.Counter()];minq=[collections.Counter(),collections.Counter()];byname={};n=0
 for one,two in itertools.zip_longest(records(a.r1),records(a.r2)):
  if one is None or two is None:raise ValueError('mate EOF disagreement')
  if one[0]!=two[0]:raise ValueError('mate name disagreement')
  n+=1
  if n>100000:raise ValueError('Bounded control inventory stops above100000 pairs')
  name=one[0];names[name]+=1;raw=one[1]+'|'+two[1];fam=hashlib.sha256(raw.encode()).hexdigest();nfam=hashlib.sha256(min(raw,rc(two[1])+'|'+rc(one[1])).encode()).hexdigest();families[fam]+=1;normalized[nfam]+=1
  byname.setdefault(name,[]).append({'full_pair_sequence_sha256':fam,'orientation_normalized_pair_sha256':nfam,'lengths':[len(one[1]),len(two[1])],'qualities_sha256':[hashlib.sha256(one[2].encode()).hexdigest(),hashlib.sha256(two[2].encode()).hexdigest()]})
  for i,(_,s,q) in enumerate([one,two]):lengths[i][len(s)]+=1;minq[i][min(ord(c)-33 for c in q)]+=1
 out={'status':'passed','paired_records':n,'unique_query_names':len(names),'repeated_query_names':{k:v for k,v in names.items() if v>1},'distinct_full_pair_sequences':len(families),'distinct_orientation_normalized_pair_sequences':len(normalized),'full_pair_family_multiplicities':dict(collections.Counter(families.values())),'name_sequence_map':byname,'length_distributions':[dict(x) for x in lengths],'whole_read_minimum_quality_distributions':[dict(x) for x in minq],'source_files':[{'path':x,'bytes':pathlib.Path(x).stat().st_size,'sha256':checksum(x)} for x in [a.r1,a.r2]],'paired_EOF_name_and_quality_checks':True,'elapsed_seconds':time.time()-t,'limitations':['Control-set sequence/name inventory only. Unique names or sequence families are not proof of independent biological molecules.','Does not itself label any sequence as a true fusion or establish mappability.']}
 pathlib.Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','paired_records','unique_query_names','distinct_full_pair_sequences','distinct_orientation_normalized_pair_sequences']}))
if __name__=='__main__':main()

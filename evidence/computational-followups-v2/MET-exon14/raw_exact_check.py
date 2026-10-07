import pathlib,json,re,gzip,collections,datetime,time
B=pathlib.Path(__file__).resolve().parent;D=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');data=json.load(open(B/'junction-check.json'))
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
lookup=collections.defaultdict(list)
for label,v in data['results'].items():
 full=v['RNA_junction_25bp_each_side']
 for arm in (15,25):
  seq=full[25-arm:25+arm]
  for orientation,p in [('+',seq),('-',rc(seq))]:lookup[p].append((label,arm,orientation))
regex=re.compile('(?=('+('|'.join(lookup))+'))')
def hit(seq,qual):
 result=set()
 for m in regex.finditer(seq):
  p=m.group(1);q=qual[m.start():m.start()+len(p)]
  if min(q,default='!')>='5':result.update((a,b) for a,b,c in lookup[p])
 return result
# Exact, reverse and low-quality controls; separate straightforward .find for each pattern.
for p,labels in lookup.items():
 expected={(a,b) for a,b,c in labels}
 assert expected<=hit('ACTGA'+p+'TTGCA','I'*(len(p)+10))
 assert not hit(p,'!'*len(p))
 q='AC'+p+p+'TA';got=hit(q,'I'*len(q));manual=set()
 for pattern,ls in lookup.items():
  if pattern in q:manual.update((a,b) for a,b,c in ls)
 assert got==manual
counts=collections.Counter();names=collections.defaultdict(list);n=0;t=time.time()
with gzip.open(D/'RNA_TN26-279853_S25.R1.fastq.gz','rt') as f,gzip.open(D/'RNA_TN26-279853_S25.R2.fastq.gz','rt') as g:
 while True:
  a=f.readline();b=g.readline()
  if not a:assert not b;break
  sa=f.readline().strip();sb=g.readline().strip();pa=f.readline();pb=g.readline();qa=f.readline().strip();qb=g.readline().strip()
  assert len(sa)==len(qa) and len(sb)==len(qb)
  found=hit(sa,qa)|hit(sb,qb)
  for k in found:counts[k]+=1;names[k].append(a.strip().split()[0])
  n+=1
assert n==23209264
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_pairs':n,'elapsed_seconds':round(time.time()-t,2),'method':'Alignment-independent exact junction pattern in either mate/orientation; quality>=20 at every marker base. Count each pair once per target/arm length. Exact30/50nt markers can miss mismatched/shorter junction reads; not clinical assay sensitivity.','public_synthetic_controls':'exact/reverse/lowquality/overlapping-regex versusstr.find checks pass','counts':{label:{str(arm):counts[(label,arm)] for arm in (15,25)} for label in data['results']},'names':{label:{str(arm):names[(label,arm)] for arm in (15,25)} for label in data['results']}}
(B/'raw-exact-check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='names'},indent=2))

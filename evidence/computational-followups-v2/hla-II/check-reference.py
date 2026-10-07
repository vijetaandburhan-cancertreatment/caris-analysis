"""Exact reference translation check; not a normal-tissue or cross-reactivity assay."""
import gzip,re,pathlib,json,hashlib,time,resource,collections
P=pathlib.Path(__file__).resolve().parent;d=json.loads((P/'patient-peptide-manifest.json').read_text());ref=P.parents[1]/'oct1-deep/neoantigen/gencode50-sensitivity/gencode.v50.pc_translations.fa.gz'
start=time.time();longs={w['mutant'] for w in d['windows']};core_meta=collections.defaultdict(lambda:{'genes':set(),'mutation_containing':False,'shared_WT':False})
for w in d['windows']:
 for i in range(len(w['mutant'])-8):
  s=w['mutant'][i:i+9];core_meta[s]['genes'].add(w['gene']);core_meta[s]['mutation_containing']|= any(i+1<=p<=i+9 for p in w['altered_positions_in_peptide1']);core_meta[s]['shared_WT']|= s in w['wt']
queries=longs|set(core_meta);prefix=collections.defaultdict(list)
for q in queries:prefix[q[:9]].append(q)
pat=re.compile('(?=('+'|'.join(sorted(prefix))+'))')
res={q:{'reference_occurrences':0,'reference_records':0,'example_records':[]} for q in queries};n=0

def scan(name,seq):
 seen=set()
 for m in pat.finditer(seq):
  for q in prefix[m.group(1)]:
   if seq.startswith(q,m.start()):
    r=res[q];r['reference_occurrences']+=1;seen.add(q)
    if len(r['example_records'])<10:r['example_records'].append({'header':name,'amino_acid_start1':m.start()+1})
 for q in seen:res[q]['reference_records']+=1
name=None;buf=[]
with gzip.open(ref,'rt') as f:
 for line in f:
  if line.startswith('>'):
   if name is not None:scan(name,''.join(buf));n+=1
   name=line.strip()[1:];buf=[]
  else:buf.append(line.strip())
 if name is not None:scan(name,''.join(buf));n+=1
for q,m in core_meta.items():m['genes']=sorted(m['genes'])
out={'reference':str(ref),'reference_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'reference_records_scanned':n,'long_peptide_count':len(longs),'core_count':len(core_meta),'novel_sequence_core_count':sum(m['mutation_containing'] for m in core_meta.values()),'long_peptide_exact_matches':{q:res[q] for q in sorted(longs) if res[q]['reference_occurrences']},'cores':{q:{**m,**res[q]} for q,m in sorted(core_meta.items())},'elapsed_seconds':time.time()-start,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'limitations':'Reference translations, not individual germline or expressed normal proteome. Nine-AA windows are possible cores, not inferred binding registers. Exact absence is not proof of tumor specificity, presentation, T-cell recognition or safety.'}
(P/'normal-reference-exact-check.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='cores'},indent=2))

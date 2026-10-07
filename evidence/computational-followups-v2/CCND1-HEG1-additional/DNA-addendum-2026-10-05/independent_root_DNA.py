from pathlib import Path
import pysam, json, hashlib, collections, datetime

W=Path(__file__).resolve().parent
L=Path.home()/'.local/share/codex/caris-analysis'
src=W/'bounded-DNA-records.bam'
truth=json.loads((W/'dna-audit.json').read_text())
f=pysam.FastaFile(str(L/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'))
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
model=rc(f.fetch('chr11',69651216,69651466).upper())+rc(f.fetch('chr3',125013072,125013572).upper())
markers={}
for offset in range(233,251):
 for n in (50,100):
  s=model[offset-n//2:offset+n//2]
  for t in (s,rc(s)): markers[t]=n
assert len(markers)==72
assert set(markers)=={x['sequence'] for x in truth['marker_definitions']}

def scan(seq,qual):
 result=set()
 for n in (50,100):
  for i in range(len(seq)-n+1):
   if seq[i:i+n] in markers:
    for q in (20,30):
     if min(qual[i:i+n])>=q: result.add((n,q))
 return result

control_count=0
for seq,n in markers.items():
 for q in (19,20,29,30):
  hits=scan('ACG'+seq+'TCA',[40]*3+[q]*n+[40]*3)
  assert ((n,20) in hits)==(q>=20)
  assert ((n,30) in hits)==(q>=30)
  control_count+=1

sets=collections.defaultdict(set); positions=collections.defaultdict(set); anchors=collections.defaultdict(set)
N=0; strictN=0; names=set(); noBQexact=set(); BAM=pysam.AlignmentFile(src,'rb')
for a in BAM:
 N+=1
 key=(a.query_name,a.get_tag('RG') if a.has_tag('RG') else '')
 names.add(key)
 seq=a.query_sequence or ''; qualities=a.query_qualities
 strict=not(a.flag & (4|256|512|1024|2048)) and a.mapping_quality>=20
 strictN+=strict
 for n in (50,100):
  if any(seq[i:i+n] in markers for i in range(len(seq)-n+1)): noBQexact.add(key)
 if qualities is None: continue
 for n,q in scan(seq,qualities):
  sets[f'all_flags_{n}nt_BQ{q}'].add(key)
  if strict:sets[f'strict_{n}nt_BQ{q}'].add(key)
 if not strict:continue
 gene='CCND1' if a.reference_name=='chr11' else 'HEG1' if a.reference_name=='chr3' else None
 if not gene:continue
 mp={};rp=a.reference_start;qp=0
 for op,k in a.cigartuples or []:
  if op in (0,7,8):
   for i in range(k):mp[rp+i]=qp+i
   rp+=k;qp+=k
  elif op in (2,3):rp+=k
  elif op in (1,4):qp+=k
  elif op not in (5,6):raise ValueError(op)
 g=truth['regions'][gene]
 for p in g['positions1']:
  q=mp.get(p-1)
  if q is not None and qualities[q]>=20:positions[(gene,str(p))].add(key)
  for n in (25,50):
   st=p-1 if gene=='CCND1' else p-n
   ix=[mp.get(k) for k in range(st,st+n)]
   if any(k is None for k in ix):continue
   if ix!=list(range(ix[0],ix[0]+n)) or min(qualities[k] for k in ix)<20:continue
   anchors[(gene,f'{p}_{n}nt','aligned_contiguous_BQ20_fragments')].add(key)
   if seq[ix[0]:ix[-1]+1]==f.fetch(g['chrom'],st,st+n).upper():anchors[(gene,f'{p}_{n}nt','reference_exact_BQ20_fragments')].add(key)
counts={k:len(sets[k]) for k in truth['exact_marker_counts']}
assert counts==truth['exact_marker_counts']
assert N==truth['final_distinct_alignment_records_including_mates']
assert strictN==truth['strict_alignment_records']
assert len(names)==truth['final_distinct_name_RG_fragments']
coverage_checked=0
for gene,cov in truth['coverage'].items():
 for p,x in cov['bases_by_equivalent_breakpoint'].items():
  assert len(positions[(gene,p)])==x['fragments_aligned_BQ20'];coverage_checked+=1
 for p,x in cov['retained_anchors'].items():
  for field in ('aligned_contiguous_BQ20_fragments','reference_exact_BQ20_fragments'):
   assert len(anchors[(gene,p,field)])==x[field];coverage_checked+=1
out={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'independent_method':'Rebuilt joined sequence from genomic FASTA; read-substring lookup instead of marker-string search; separate CIGAR coordinate traversal instead of get_aligned_pairs; all bounded BAM records checked. Original whole DNA was not rescanned.','records':N,'strict_records':strictN,'name_RG_fragments':len(names),'exact_marker_counts':counts,'any_exact_marker_ignoring_BQ_fragments':len(noBQexact),'synthetic_quality_threshold_cases':control_count,'coverage_fields_reproduced':coverage_checked,'marker_set_rebuilt_from_genomic_FASTA_matches':True,'bounded_BAM_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scope':'Negative only within saved locus/mate records and these RNA-derived exact markers; not genomic rearrangement exclusion.'}
(W/'independent-root-DNA.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))

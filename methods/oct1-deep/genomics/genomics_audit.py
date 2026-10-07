"""Bounded research audit: no clinical calls, no original edits, no sequence export."""
import pathlib,json,csv,collections,statistics,hashlib,time,datetime
import pysam
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent
SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
OLD=ROOT/'work/oct1-analysis'
GENES={'BAP1','CDKN2A','CDKN2B','MTAP','NF2','SETD2','TP53','TERT','LATS1','LATS2','PBRM1','TRAF7','SUFU','PTEN','RB1','SMARCA4','SMARCB1','STAG2','SMARCA2','COL2A1','PRDM6'}
def dump(name,x): (OUT/name).write_text(json.dumps(x,indent=2)+'\n')
def tsv(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,keys,delimiter='\t');w.writeheader();w.writerows(rows)

# Revalidate resident inventory against recent full CRC64/SHA verification.
previous=json.loads((OLD/'resident-copies-verification.json').read_text())
reval=[]
for f in previous['files']:
 p=SRC/f['name'];st=p.stat();r={'name':p.name,'bytes_now':st.st_size,'bytes_match_previous_verified':st.st_size==f['size_bytes'],'previous_status':f['status'],'previous_sha256':f['calculated']['sha256'],'mtime_utc':datetime.datetime.fromtimestamp(st.st_mtime,datetime.timezone.utc).isoformat()}
 if st.st_size<2000000:
  h=hashlib.sha256(p.read_bytes()).hexdigest();r.update(sha256_now=h,sha256_matches_previous_verified=h==f['calculated']['sha256'])
 reval.append(r)
assert len(reval)==9 and all(r['bytes_match_previous_verified'] and r['previous_status']=='verified' and r.get('sha256_matches_previous_verified',True) for r in reval)
dump('resident-revalidation.json',{'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'All nine existence/size vs recently full-source-verified resident manifest; fresh full SHA256 of three small files. Large files not rehashed; BAM indexed reads and completed FASTQ scans are separate checks.','files':reval})

v=[]
for n,line in enumerate((SRC/'DNA_TN26-279853.vcf').open(),1):
 if line.startswith('#'):continue
 a=line.rstrip().split('\t');info=dict(z.split('=',1) if '=' in z else (z,True) for z in a[7].split(';'));sf=dict(zip(a[8].split(':'),a[9].split(':')))
 v.append({'line':n,'chrom':a[0],'pos':int(a[1]),'id':a[2],'ref':a[3],'alt':a[4],'filter':a[6],'gene':info.get('GI'),'consequence':info.get('FC'),'protein':info.get('PC'),'coding':info.get('DC'),'clinical':info.get('CI'),'DP':int(info['DP']),'VAF':float(sf['VF']),'AD':sf['AD'],'SA':sf['SA'],'GT':sf['GT']})
assert len(v)==2843
selected=list(csv.DictReader((ROOT/'outputs/caris-analysis/Caris-candidate-evidence.tsv').open(),delimiter='\t'))
for s in selected:
 hits=[r for r in v if (r['chrom'],str(r['pos']),r['ref'],r['alt'])==(s['hg38_chrom'],s['position_1based'],s['ref'],s['alt'])]
 assert len(hits)==1 and hits[0]['filter']=='.' and abs(hits[0]['VAF']-float(s['Caris_VCF_tumor_vaf']))<1e-9
dump('candidate-revalidation.json',{'selected_count':len(selected),'all21_source_loci_filter_vaf_exact_match':True,'VCF_records':len(v)})
tsv('target-gene-vcf-records.tsv',[r for r in v if r['gene'] in GENES])
x=json.loads((OLD/'variants.xlsx-records.json').read_text())
cols=['source_row','Test','Technology','Biomarker','Test Result','CNA Value','NGS Transcript','NGS Protein Change','NGS Interpretation','Fusion Gene 1','Fusion Gene 2','Breakpoint','Total Uniq Read Pairs']
tsv('target-gene-workbook-tests.tsv',[{k:r.get(k,'') for k in cols} for r in x if r['Biomarker'] in GENES])
tsv('indeterminate-tests.tsv',[{k:r.get(k,'') for k in cols} for r in x if 'indeter' in r['Test Result'].lower()])

# Descriptive germline-like marker imbalance. dbSNP is not proof of germline.
# Never promote Benign/rs markers into cancer driver candidates.
markers=[r for r in v if r['id'].startswith('rs') and len(r['ref'])==len(r['alt'])==1 and .1<r['VAF']<.9 and r['DP']>=100 and set(r['filter'].split(';')) <= {'.','Benign','rs'}]
regions={'chr3_0to80Mb':('chr3',0,80000000),'chr3_100Mbtoend':('chr3',100000000,200000000),'chr9_all':('chr9',0,140000000),'chr17_0to20Mb':('chr17',0,20000000),'chr17_30Mbtoend':('chr17',30000000,85000000)}
summary={}
for name,(ch,st,en) in regions.items():
 rows=[r for r in markers if r['chrom']==ch and st<=r['pos']<en]
 for r in rows:r.setdefault('region',name)
 m=[min(r['VAF'],1-r['VAF']) for r in rows]
 summary[name]={'n':len(rows),'median_minor_allele_fraction':statistics.median(m) if m else None,'range_minor_allele_fraction':[min(m),max(m)] if m else [],'note':'Selected dbSNP biallelic markers only; sparse gene-targeted and ascertainment-biased. Neither normal heterozygosity nor allele-specific copy number established.'}
dump('allelic-imbalance-descriptive-summary.json',summary)
chosen=[r for r in markers if 'region' in r]
tsv('allelic-imbalance-markers.tsv',chosen)

def count_snv(bam,ch,p):
 reads=collections.Counter();strand=collections.Counter();frags=collections.defaultdict(set);clean=collections.defaultdict(set)
 for r in bam.fetch(ch,p-1,p):
  if r.flag & (4|256|512|1024|2048) or r.mapping_quality<20:continue
  q=next((q for q,rp in r.get_aligned_pairs(matches_only=True) if rp==p-1),None)
  if q is None or r.query_qualities is None or r.query_qualities[q]<20:continue
  b=r.query_sequence[q];reads[b]+=1;strand[b+('_reverse' if r.is_reverse else '_forward')]+=1
  key=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name);frags[key].add(b)
  if not any(op==4 for op,n in r.cigartuples) and min(q-r.query_alignment_start,r.query_alignment_end-1-q)>=5:clean[b].add(key)
 fc=collections.Counter(next(iter(a)) if len(a)==1 else 'discordant' for a in frags.values())
 return {'read_counts':dict(reads),'fragment_counts':dict(fc),'read_strands':dict(strand),'strict_no_clip_no_end_fragment_counts':{k:len(s) for k,s in clean.items()}}

# Bound independent DNA recount to 13 informative markers + two TERT hotspots.
recount_ids={'rs2272125','rs2228001','rs2279017','rs2227998','rs200001131','rs1385816','rs2305037','rs3218651','rs487848','rs1042522','rs1800369','rs2230722','rs2229974'}
calls=[]
with pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(OLD/'DNA_TN26-279853.bam.bai')) as bam:
 for r in markers:
  if r['id'] in recount_ids: calls.append({**r,**count_snv(bam,r['chrom'],r['pos'])})
 for p,label in [(1295113,'TERT_C228T_c.-124C>T'),(1295135,'TERT_C250T_c.-146C>T')]:
  calls.append({'chrom':'chr5','pos':p,'label':label,'ref':'G','alt':'A','source_for_coordinates':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11914184/',**count_snv(bam,'chr5',p)})
dump('bounded-DNA-allele-recounts.json',{'method':'MAPQ/BQ>=20, exclude unmapped/secondary/supplementary/QCfail/duplicates; query-name+RG mate collapse, disagreement excluded. TERT ref G > A on forward genome (coding strand C>T).','calls':calls})
print(json.dumps({'revalidation':'passed','selected_candidates':len(selected),'regional_marker_summary':summary,'bounded_recounts':len(calls)},indent=2),flush=True)

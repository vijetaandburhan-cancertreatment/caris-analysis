from pathlib import Path
import csv,json,bisect,collections,time,hashlib,gzip,pysam
R=Path(__file__).resolve().parent;P=R.parent/'oct5-population-panel-concordance-v1';O=R/'comparison';O.mkdir(exist_ok=True)
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,rs):
 with p.open('w') as f:
  if not rs:return
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def js(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
assert json.loads((R/'patient-full-pass/run.json').read_text())['status']=='COMPLETE'
cohort=rows(R/'baseline-fixed-4015.tsv');assert len(cohort)==4015
lookup={(x['chrom'],int(x['pos1'])-1):x for x in cohort};sites=collections.defaultdict(list)
for ch,pos in lookup:sites[ch].append(pos)
for ch in sites:sites[ch].sort()
groups={k:{} for k in lookup};clean={k:{} for k in lookup};reads={k:collections.Counter() for k in lookup};mask={b:1<<i for i,b in enumerate('ACGT')};stats=collections.Counter();t=time.time();exceptions={(x['chrom'],int(x['pos1'])-1) for x in rows(R/'frozen-exception-sites.tsv')}
bam=R/'patient-full-pass/cohort-and-exception-records.bam'
with pysam.AlignmentFile(str(bam),'rb') as f:
 hdr=f.header.to_dict();refs=f.references
 for rd in f.fetch(until_eof=True):
  stats['captured_records_seen']+=1
  if rd.flag&0xF0C or not rd.flag&2 or rd.mapping_quality<30:continue
  ch=refs[rd.reference_id];ps=sites.get(ch,());idx=bisect.bisect_left(ps,rd.reference_start)
  if idx==len(ps) or ps[idx]>=rd.reference_end:continue
  seq=rd.query_sequence;quals=rd.query_qualities
  if seq is None or quals is None:continue
  rp=rd.reference_start;qp=0;no_soft=not any(op==4 for op,n in rd.cigartuples)
  for op,n in rd.cigartuples:
   if op in (0,7,8):
    lo=bisect.bisect_left(ps,rp);hi=bisect.bisect_left(ps,rp+n)
    for pos in ps[lo:hi]:
     q=qp+pos-rp
     if quals[q]<25 or seq[q].upper() not in mask:continue
     k=(ch,pos);b=seq[q].upper();g=groups[k];g[rd.query_name]=g.get(rd.query_name,0)|mask[b];reads[k][b]+=1
     if no_soft and q-rd.query_alignment_start>=5 and rd.query_alignment_end-1-q>=5:
      g=clean[k];g[rd.query_name]=g.get(rd.query_name,0)|mask[b]
    rp+=n;qp+=n
   elif op in(1,4):qp+=n
   elif op in(2,3):rp+=n
   elif op in(5,6):pass
   else:raise RuntimeError(op)
result=[];transition=collections.Counter();fixed=collections.Counter();new_events=[]
for k,r0 in lookup.items():
 r={k:v for k,v in r0.items() if k in ['chrom','pos1','id','ref','alt','overlap_gene_symbols','DNA_category','DNA_ACGT_depth','DNA_fragment_ref','DNA_fragment_alt','DNA_fragment_other','RNA_ACGT_depth','RNA_fragment_ref','RNA_fragment_alt','RNA_fragment_other','RNA_same_allele_fraction','retain90','retain98','frozen_spaced50kb']}
 for label,source in [('new',groups),('new_clean',clean)]:
  ct=collections.Counter(source[k].values());ac={b:ct[m] for b,m in mask.items()};depth=sum(ac.values())
  for b in 'ACGT':r[label+'_'+b]=ac[b]
  r[label+'_depth']=depth;r[label+'_conflicting_names']=sum(v for m,v in ct.items() if m not in mask.values())
  r[label+'_ref']=ac[r['ref']];r[label+'_alt']=ac[r['alt']];r[label+'_other']=depth-ac[r['ref']]-ac[r['alt']]
 r['new_callable20']=int(r['new_depth']>=20);r['new_zero_ACGT']=int(r['new_depth']==0);r['new_same_fixed_DNA_fraction']='';r['new_retain90_with_coverage20']=0;r['new_retain98_with_coverage20']=0;r['new_opposite98_with_coverage20']=0
 r['new_read_ref']=reads[k][r['ref']];r['new_read_alt']=reads[k][r['alt']]
 fixed['all_original_sites']+=1;fixed['new_callable20']+=r['new_callable20'];fixed['new_zero_ACGT']+=r['new_zero_ACGT'];fixed['new_below20_nonzero']+=int(0<r['new_depth']<20)
 cat=r['DNA_category']
 if cat in ('ref_dominant','alt_dominant'):
  fixed['original_DNA_dominant']+=1;b=r['ref'] if cat=='ref_dominant' else r['alt'];other=r['alt'] if cat=='ref_dominant' else r['ref'];d=r['new_depth'];same=r['new_'+b]
  r['new_same_fixed_DNA_fraction']=same/d if d else ''
  r['new_retain90_with_coverage20']=int(d>=20 and same*10>=d*9);r['new_retain98_with_coverage20']=int(d>=20 and same*100>=d*98);r['new_opposite98_with_coverage20']=int(d>=20 and r['new_'+other]*100>=d*98)
  fixed['dominant_new_callable20']+=int(d>=20);fixed['dominant_new_zero']+=int(d==0);fixed['dominant_new_below20_nonzero']+=int(0<d<20);fixed['dominant_same90_fixed_denominator2873']+=r['new_retain90_with_coverage20'];fixed['dominant_same98_fixed_denominator2873']+=r['new_retain98_with_coverage20'];fixed['dominant_opposite98']+=r['new_opposite98_with_coverage20']
  oldstate='retain90' if int(r['retain90']) else 'below90';newstate='below20' if d<20 else 'retain90' if r['new_retain90_with_coverage20'] else 'below90';transition[oldstate+'__'+newstate]+=1
 result.append(r)
write(O/'fixed4015-comparison.tsv',result);write(O/'all16-exception-comparison.tsv',[r for r in result if (r['chrom'],int(r['pos1'])-1) in exceptions]);write(O/'new_coverage_losses.tsv',[r for r in result if not r['new_callable20']])
# Exact original DNA per-base evidence was independently produced in the previous milestone.
dna={(r['chrom'],int(r['pos1'])-1):r for r in rows(P/'independent-count-audit/DNA-direct-counts.tsv')}
for r in result:
 k=(r['chrom'],int(r['pos1'])-1);dn=dna[k];dd=sum(int(dn['baseline_'+b]) for b in 'ACGT');nd=r['new_depth']
 if min(dd,nd)<50:continue
 for b in 'ACGT':
  old=int(dn['baseline_'+b]);new=r['new_'+b]
  if old*100<=2*dd and new*10>=nd and new>=10:new_events.append({'chrom':k[0],'pos1':k[1]+1,'gene':r['overlap_gene_symbols'],'base':b,'DNA_count':old,'DNA_depth':dd,'new_RNA_count':new,'new_RNA_depth':nd,'is_frozen16_exception':k in exceptions})
write(O/'new_strong_allele_flags.tsv',new_events)
with gzip.open(O/'new-exception-QNAME-alleles.json.gz','wt') as f:json.dump({f'{ch}:{pos+1}':{'baseline_filter':groups[(ch,pos)],'clean_filter':clean[(ch,pos)]} for ch,pos in sorted(exceptions)},f)
summary={'status':'PRIMARY_COMPARISON_COMPLETE_PENDING_INDEPENDENT_AUDIT','fixed_cohort':dict(fixed),'dominant_transition':dict(transition),'new_strong_events':new_events,'new_strong_loci':len({(r['chrom'],r['pos1']) for r in new_events}),'baseline_denominators_preserved':True,'original_DNA_retained_unchanged':True,'elapsed_seconds':time.time()-t,'captured_record_stats':dict(stats),'new_header':hdr,'captured_BAM_sha256':sha(bam),'cohort_sha256':sha(R/'baseline-fixed-4015.tsv'),'limits':'Matched count filters; different upstream mapping/reference/parameters. Coverage loss remains in fixed denominators. Capture holds onlySNP-overlap records for ordinary names and everyemittedrecord for frozen3426names. D/N-only bookkeeping notcomplete; noformal clinicalperformance claim.'};js(O/'summary.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['new_header','new_strong_events']},indent=2))

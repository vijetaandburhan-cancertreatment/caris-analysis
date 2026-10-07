from pathlib import Path
import json,csv,gzip,pysam,collections,hashlib
R=Path(__file__).resolve().parent;O=R/'name-comparison';O.mkdir(exist_ok=True)
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def save(p,data):p.write_text(json.dumps(data,indent=2)+'\n')
def table(p,rs):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
base_bits={b:1<<i for i,b in enumerate('ACGT')};single={v:k for k,v in base_bits.items()};sites=rows(R/'frozen-exception-sites.tsv');site_by_chr=collections.defaultdict(list)
for s in sites:site_by_chr[s['chrom']].append((int(s['pos1'])-1,s))
names=set((R/'frozen-exception-qnames.txt').read_text().splitlines())
assert json.loads((R/'patient-full-pass/run.json').read_text())['status']=='COMPLETE'
assert json.loads((R/'original-complete-names/summary.json').read_text())['status']=='COMPLETE'
def inspect(path):
 records=collections.defaultdict(list);calls=collections.defaultdict(dict);overlap=collections.defaultdict(set);gaps=collections.defaultdict(set);observations=collections.defaultdict(list)
 with pysam.AlignmentFile(str(path),'rb') as f:
  for rd in f.fetch(until_eof=True):
   if rd.query_name not in names:continue
   records[rd.query_name].append(rd)
   if rd.is_unmapped or not rd.query_sequence:continue
   options=[(p,s) for p,s in site_by_chr.get(rd.reference_name,[]) if rd.reference_start<=p<rd.reference_end]
   if not options:continue
   amap={rp:qp for qp,rp in rd.get_aligned_pairs(matches_only=False) if rp is not None}
   for p,s in options:
    k=f'{s["chrom"]}:{p+1}';q=amap.get(p)
    if q is None:gaps[k].add(rd.query_name);continue
    overlap[k].add(rd.query_name)
    if rd.flag&0xF0C or not rd.flag&2 or rd.mapping_quality<30 or rd.query_qualities is None or rd.query_qualities[q]<25:continue
    b=rd.query_sequence[q].upper()
    if b not in base_bits:continue
    calls[k][rd.query_name]=calls[k].get(rd.query_name,0)|base_bits[b]
    physical=len(rd.query_sequence)-1-q if rd.is_reverse else q
    observations[(k,rd.query_name)].append({'base':b,'read1':rd.is_read1,'read2':rd.is_read2,'forward_query_offset':physical,'query_length':len(rd.query_sequence),'has_hard_clip':any(op==5 for op,n in rd.cigartuples),'forward_sequence_sha256':hashlib.sha256(rd.get_forward_sequence().encode()).hexdigest(),'flag':rd.flag,'position1':rd.reference_start+1,'cigar':rd.cigarstring})
 return records,calls,overlap,gaps,observations
oldr,oldc,oldo,oldg,oldobs=inspect(R/'original-complete-names/original-3426-names.bam')
newr,newc,newo,newg,newobs=inspect(R/'patient-full-pass/cohort-and-exception-records.bam')
assert set(oldr)==names
with gzip.open(R/'comparison/new-exception-QNAME-alleles.json.gz','rt') as f:allnew=json.load(f)
result=[];transitions=[];record_summaries=[];coordinate_transitions=[]
for s in sites:
 k=f'{s["chrom"]}:{s["pos1"]}';counts=collections.Counter(oldc[k].values());observed=','.join(str(counts[base_bits[b]]) for b in 'ACGT');assert observed==s['baseline_RNA_ACGT'],(k,observed,s['baseline_RNA_ACGT'])
 assert {n:m for n,m in allnew[k]['baseline_filter'].items() if n in names}==newc[k]
 if s['DNA_category']=='alt_dominant':expected=s['alt']
 elif s['DNA_category']=='ref_dominant':expected=s['ref']
 else:expected=None
 unexpected_bases=[s['unexpected_base']] if s['review_class']=='MAPQ60_only_not_baseline_flag' else [b for b in 'ACGT' if b!=expected and counts[base_bits[b]]]
 for b in unexpected_bases:
  supporters={n for n,m in oldc[k].items() if m==base_bits[b]};outcomes=collections.Counter()
  for n in sorted(supporters):
   m=newc[k].get(n,0);rr=newr.get(n,[])
   if m==base_bits[b]:status='retains_same_unexpected_base'
   elif m in single:status='now_'+single[m]+('_DNA_expected' if single[m]==expected else '_other_base')
   elif m:status='conflicting_qualified_bases'
   elif not rr:status='no_emitted_record'
   elif n in newo[k]:status='site_base_overlap_but_fails_matched_filters'
   elif n in newg[k]:status='site_reference_span_but_no_query_base'
   elif any(not r.is_unmapped for r in rr):status='mapped_without_base_at_original_site'
   else:status='only_unmapped_emitted_records'
   outcomes[status]+=1
   transitions.append({'site':k,'gene':s['overlap_gene_symbols'],'review_class':s['review_class'],'original_unexpected_base':b,'QNAME':n,'status':status,'new_qualified_mask':m,'old_record_count':len(oldr[n]),'new_record_count':len(rr),'new_primary_mapped_records':sum(not r.is_unmapped and not r.is_secondary and not r.is_supplementary for r in rr),'new_secondary_records':sum(r.is_secondary for r in rr),'new_supplementary_records':sum(r.is_supplementary for r in rr),'new_unmapped_records':sum(r.is_unmapped for r in rr),'new_locations':';'.join(sorted({f'{r.reference_name}:{r.reference_start+1}:{r.cigarstring}:flag{r.flag}:MQ{r.mapping_quality}' for r in rr if not r.is_unmapped}))})
   # Trace the same physical query offset from qualifying old unexpected observations.
   old_obs={(x['read1'],x['read2'],x['forward_query_offset'],x['query_length'],x['has_hard_clip'],x['forward_sequence_sha256']) for x in oldobs[(k,n)] if x['base']==b}
   for read1,read2,physical,qlen,old_hardclip,old_seq_hash in sorted(old_obs):
    for r in rr:
     if (r.is_read1,r.is_read2)!=(read1,read2) or r.query_sequence is None:continue
     if old_hardclip or any(op==5 for op,n in (r.cigartuples or [])) or len(r.query_sequence)!=qlen or hashlib.sha256(r.get_forward_sequence().encode()).hexdigest()!=old_seq_hash:
      coordinate_transitions.append({'site':k,'QNAME':n,'old_base':b,'mate':'R1' if read1 else 'R2','forward_query_offset':physical,'new_flag':r.flag,'new_MAPQ':r.mapping_quality,'new_coordinate':'not_traced_hardclip_length_or_forward_sequence_difference','new_CIGAR':r.cigarstring or '*'});continue
     q=qlen-1-physical if r.is_reverse else physical
     amap=dict(r.get_aligned_pairs(matches_only=False)) if not r.is_unmapped else {}
     rp=amap.get(q);coordinate='unmapped_or_unaligned_query_base' if rp is None else f'{r.reference_name}:{rp+1}'
     coordinate_transitions.append({'site':k,'QNAME':n,'old_base':b,'mate':'R1' if read1 else 'R2','forward_query_offset':physical,'new_flag':r.flag,'new_MAPQ':r.mapping_quality,'new_coordinate':coordinate,'new_CIGAR':r.cigarstring or '*'})
  result.append({'site':k,'gene':s['overlap_gene_symbols'],'review_class':s['review_class'],'original_unexpected_base':b,'original_supporters':len(supporters),'new_all_names_same_base':sum(m==base_bits[b] for m in allnew[k]['baseline_filter'].values()),'old_supporter_outcomes':dict(outcomes)})
for n in sorted(names):
 for label,rs in [('original',oldr.get(n,[])),('new',newr.get(n,[]))]:
  record_summaries.append({'QNAME':n,'assay_alignment':label,'records':len(rs),'read1_records':sum(r.is_read1 for r in rs),'read2_records':sum(r.is_read2 for r in rs),'unmapped_records':sum(r.is_unmapped for r in rs),'secondary_records':sum(r.is_secondary for r in rs),'supplementary_records':sum(r.is_supplementary for r in rs),'locations':';'.join(sorted({f'{r.reference_name}:{r.reference_start+1}:{r.cigarstring}:flag{r.flag}:MQ{r.mapping_quality}' for r in rs if not r.is_unmapped}))})
a='chr8:100705591';bb='chr8:100705604';og={n for n,m in oldc[a].items() if m==base_bits['G']};oc={n for n,m in oldc[bb].items() if m==base_bits['C']};ng={n for n,m in allnew[a]['baseline_filter'].items() if m==base_bits['G']};nc={n for n,m in allnew[bb]['baseline_filter'].items() if m==base_bits['C']}
pab={'original_591_G_names':len(og),'original_604_C_names':len(oc),'original_shared':len(og&oc),'new_591_G_names_all':len(ng),'new_604_C_names_all':len(nc),'new_shared_all':len(ng&nc),'old_shared_still_both':len(og&oc&ng&nc),'old_shared_retained_591_G':len(og&oc&ng),'old_shared_retained_604_C':len(og&oc&nc),'new_shared_not_in_old_shared':len((ng&nc)-(og&oc))}
summary={'status':'PRIMARY_NAME_COMPARISON_COMPLETE_PENDING_INDEPENDENT_AUDIT','original_names':len(oldr),'new_frozen_names_seen':len(newr),'missing_new_frozen_names':sorted(names-set(newr)),'all16_original_counts_match_frozen':True,'new_frozen_name_masks_match_primary_full_counts':True,'unexpected_supporter_transitions':result,'PABPC1_linkage':pab,'limitations':'Same raw reads, changed full alignment pipeline. New placement is not necessarily correct. All emitted records are preserved for frozen names; not all hypothetical locations. New names outside frozen3426 have only cohort-overlap records. Auxiliary physical-offset tracing requires no hard clips in either record and exact forward-oriented query sequence equality; other cases are explicitly not traced. Shared QNAME support alone is not proof of two alleles on one record or independent molecules.'}
save(O/'summary.json',summary);table(O/'unexpected-name-transitions.tsv',transitions);table(O/'all3426-name-record-inventory.tsv',record_summaries);table(O/'physical-query-position-transitions.tsv',coordinate_transitions)
print(json.dumps(summary,indent=2))

#!/usr/bin/env python3
"""Compare exact read names and all retained alignments, including filtered ones."""
import pathlib,json,collections,hashlib,csv,pysam
W=pathlib.Path(__file__).resolve().parent;C=pathlib.Path.home()/'.local/share/codex/caris-analysis'
names=json.loads((W/'RASA1-difference-names.json').read_text());wanted=set(names['lost_alternate']+names['gained_alternate'])
s=json.loads((W/'new-D8/summary.json').read_text());ctx=s['variant_contexts']['RASA1'];g=s['gene_intervals']['RASA1']
left,right=ctx['start0']-10,ctx['end0']+10
def describe(r):
 x=r.reference_start;q=0;mapping={}
 for op,n in r.cigartuples or []:
  if op in [0,7,8]:
   for k in range(n):mapping[x+k]=q+k
  if op in [0,2,3,7,8]:x+=n
  if op in [0,1,4,7,8]:q+=n
 causes=[]
 for bit,name in [(4,'unmapped'),(256,'secondary'),(512,'QCfail'),(1024,'duplicate'),(2048,'supplementary')]:
  if r.flag&bit:causes.append(name)
 if r.mapping_quality<20:causes.append('MAPQ<20')
 if r.has_tag('NH') and r.get_tag('NH')!=1:causes.append('NH!=1')
 a,z=mapping.get(left),mapping.get(right-1);window=None;minq=None;allele=None
 if a is None or z is None:coverage='missing_left_anchor' if a is None and z is not None else 'missing_right_anchor' if z is None and a is not None else 'missing_both_anchors'
 elif z<a:coverage='nonincreasing_query_anchors'
 else:
  coverage='both_anchors';window=r.query_sequence[a:z+1];minq=min(r.query_qualities[a:z+1]) if r.query_qualities is not None else None
  allele='alternate' if window==ctx['alternate'] else 'reference' if window==ctx['reference'] else 'other'
  if minq is None or minq<20:coverage='low_or_missing_BQ'
 return {'mate':1 if r.is_read1 else 2 if r.is_read2 else 0,'flag':r.flag,'chrom':r.reference_name,'start0':r.reference_start,'end0':r.reference_end,'cigar':r.cigarstring,'MAPQ':r.mapping_quality,'NH':r.get_tag('NH') if r.has_tag('NH') else None,'alignment_filter_reasons':causes,'window_coverage':coverage,'window_allele_before_filter':allele,'window_min_BQ':minq,'window_sequence':window,'read_length':r.query_length,'forward_read_sha256':hashlib.sha256(r.get_forward_sequence().encode()).hexdigest(),'forward_quality_sha256':hashlib.sha256(bytes(r.get_forward_qualities())).hexdigest() if r.query_qualities is not None else None,'passes_strict_alternate':not causes and coverage=='both_anchors' and allele=='alternate'}
data={name:{'old':[],'new':[]} for name in wanted}
oldbam=pysam.AlignmentFile(str(C/'TN26-279853/RNA_TN26-279853.bam'),'rb',index_filename=str(W.parents[2]/'oct1-analysis/RNA_TN26-279853.bam.bai'))
for r in oldbam.fetch(g['chrom'],g['start0'],g['end0']):
 if r.query_name in wanted:data[r.query_name]['old'].append(describe(r))
newbam=pysam.AlignmentFile(s['bam'],'rb')
for r in newbam.fetch(until_eof=True):
 if r.query_name in wanted:data[r.query_name]['new'].append(describe(r))
rows=[];counts=collections.Counter()
for name,d in sorted(data.items()):
 lost=name in names['lost_alternate'];direction='lost' if lost else 'gained'
 oldalt=[r for r in d['old'] if r['passes_strict_alternate']];newalt=[r for r in d['new'] if r['passes_strict_alternate']]
 assert bool(oldalt)==lost and bool(newalt)!=lost
 newseqalt=[r for r in d['new'] if r['window_allele_before_filter']=='alternate']
 if lost:
  if not d['new']:reason='No new alignment retained in four target intervals'
  elif newseqalt and all(r['alignment_filter_reasons'] for r in newseqalt):reason='Alternate sequence still aligns locally, but alignment uniqueness/flag filter excludes it'
  elif any(r['window_coverage']=='low_or_missing_BQ' for r in newseqalt):reason='Alternate local sequence fails window base quality'
  elif not newseqalt:reason='New retained CIGAR/placement no longer covers exact alternate window'
  else:reason='Other; inspect per-read details'
 else:reason='New alignment restores a strict alternate window absent from old qualifying alignments'
 counts[direction+': '+reason]+=1
 shared=[]
 for o in d['old']:
  for n in d['new']:
   if o['mate']==n['mate'] and o['forward_read_sha256']==n['forward_read_sha256']:
    shared.append({'mate':o['mate'],'same_forward_sequence':True,'same_forward_qualities':o['forward_quality_sha256']==n['forward_quality_sha256']})
 d.update(direction=direction,explanation=reason,same_read_sequence_checks=shared)
 rows.append({'query_name':name,'direction':direction,'explanation':reason,'old_alignment_count':len(d['old']),'new_target_alignment_count':len(d['new']),'old_strict_alt_reads':len(oldalt),'new_strict_alt_reads':len(newalt),'new_alt_sequence_reads_before_filter':len(newseqalt),'new_alt_filter_reasons':';'.join(sorted({','.join(x['alignment_filter_reasons']) for x in newseqalt})),'old_CIGARs':';'.join(r['cigar'] or '*' for r in d['old']),'new_CIGARs':';'.join(r['cigar'] or '*' for r in d['new']),'new_MAPQs':';'.join(str(r['MAPQ']) for r in d['new']),'new_NHs':';'.join(str(r['NH']) for r in d['new'])})
report={'old_alternate_fragments':75,'new_alternate_fragments':61,'lost_alternate_names':len(names['lost_alternate']),'gained_alternate_names':len(names['gained_alternate']),'net_change':-14,'classification':dict(counts),'names':data,'limits':['Same specimen and read pool: these are computational eligibility differences, not longitudinal biology.','New target BAM includes all alignments retained in four gene intervals, including low-MAPQ/secondary/supplementary reads, but cannot reveal placements elsewhere.','Counts group read names, not distinct original molecules.']}
(W/'RASA1-read-difference-details.json').write_text(json.dumps(report,indent=2)+'\n')
with (W/'RASA1-read-difference-annotated.tsv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
print(json.dumps({k:v for k,v in report.items() if k!='names'},indent=2))

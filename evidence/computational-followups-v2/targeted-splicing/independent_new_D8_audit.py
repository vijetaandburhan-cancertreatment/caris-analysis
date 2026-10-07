#!/usr/bin/env python3
"""Read-only independent CIGAR and fragment audit; does not import production code."""
import pathlib,json,csv,collections,hashlib,pysam
W=pathlib.Path(__file__).resolve().parent
S=json.loads((W/'new-D8/summary.json').read_text())
cache=pathlib.Path.home()/'.local/share/codex/caris-analysis'
b=pysam.AlignmentFile(S['bam'],'rb');all_records=list(b.fetch(until_eof=True))
fa=pysam.FastaFile(str(cache/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'))
answer={};observed_junctions={};odd_mates=[]
for gene,g in S['gene_intervals'].items():
 ctx=S['variant_contexts'].get(gene);frags={};counts=collections.Counter();flagcounts=collections.Counter();allele_N=0;cases=[]
 for r in (r for r in all_records if not r.is_unmapped and r.reference_name==g['chrom'] and r.reference_start<g['end0'] and r.reference_end>g['start0']):
  counts['all_alignments']+=1
  for bit,label in [(4,'unmapped'),(256,'secondary'),(512,'QCfail'),(1024,'duplicate'),(2048,'supplementary')]:
   if r.flag&bit:flagcounts[label]+=1
  if r.mapping_quality<20:flagcounts['MAPQ_below20']+=1
  if r.has_tag('NH') and r.get_tag('NH')!=1:flagcounts['NH_not1']+=1
  if r.flag&3844 or r.mapping_quality<20 or (r.has_tag('NH') and r.get_tag('NH')!=1):continue
  counts['qualified_alignments']+=1
  key=(r.get_tag('RG') if r.has_tag('RG') else '')+'|'+r.query_name
  d=frags.setdefault(key,{'alleles':set(),'junctions':set(),'direct':set(),'records':[]})
  d['records'].append((r.flag,r.reference_start,r.next_reference_id,r.next_reference_start,r.cigarstring,r.template_length))
  x=r.reference_start;q=0;mapping={};junctions=[];skips=[]
  for i,(op,n) in enumerate(r.cigartuples or []):
   if op in [0,7,8]:
    for k in range(n):mapping[x+k]=q+k
   if op==3:
    junctions.append((x,x+n,i,q));skips.append((x,x+n))
   if op in [0,2,3,7,8]:x+=n
   if op in [0,1,4,7,8]:q+=n
  js=set()
  for left,right,i,jq in junctions:
   # Independent contiguous query/reference anchor reconstruction, with the same
   # immediate match-operation requirement as production.
   if not(0<i<len(r.cigartuples)-1):continue
   before,after=r.cigartuples[i-1],r.cigartuples[i+1]
   if before[0] not in [0,7,8] or after[0] not in [0,7,8] or min(before[1],after[1])<12:continue
   if left-12<g['start0'] or right+12>g['end0'] or r.query_qualities is None:continue
   coords=list(range(left-12,left))+list(range(right,right+12))
   qq=[mapping.get(pos) for pos in coords]
   if qq!=list(range(jq-12,jq+12)):continue
   if min(r.query_qualities[z] for z in qq)<20:continue
   obs=''.join(r.query_sequence[z] for z in qq)
   expected=fa.fetch(g['chrom'],left-12,left).upper()+fa.fetch(g['chrom'],right,right+12).upper()
   if obs==expected:js.add((left,right))
  d['junctions'].update(js)
  if ctx:
   left,right=ctx['start0']-10,ctx['end0']+10
   a,z=mapping.get(left),mapping.get(right-1)
   if a is not None and z is not None and z>=a and r.query_qualities is not None and min(r.query_qualities[a:z+1])>=20:
    seq=r.query_sequence[a:z+1];state='alternate' if seq==ctx['alternate'] else 'reference' if seq==ctx['reference'] else 'other'
    d['alleles'].add(state)
    if state=='alternate':
     d['direct'].update(js)
     if any(l<right and u>left for l,u in skips):allele_N+=1
 counts_by_state=collections.Counter();jmap=collections.defaultdict(lambda:{'all':set(),'alt':set(),'direct':set()})
 for key,d in frags.items():
  state=next(iter(d['alleles'])) if len(d['alleles'])==1 else 'conflict' if d['alleles'] else 'uncovered';counts_by_state[state]+=1
  records=d['records']
  if len(records)>2 or (len(records)==2 and bool(records[0][0]&64)==bool(records[1][0]&64)):
   odd_mates.append({'gene':gene,'key':key,'records':records})
  if len(records)==2:
   a,z=records
   if a[3]!=z[1] or z[3]!=a[1]:counts['two_primary_mate_coordinate_mismatches']+=1
  for j in d['junctions']:
   jmap[j]['all'].add(key)
   if state=='alternate':jmap[j]['alt'].add(key)
  if state=='alternate':
   for j in d['direct']:jmap[j]['direct'].add(key)
 observed_junctions.update({(gene,*j):{k:len(v) for k,v in sets.items()} for j,sets in jmap.items()})
 answer[gene]={'alignment_counts':dict(counts),'excluded_flag_counts_nonexclusive':dict(flagcounts),'fragment_allele_states':dict(counts_by_state),'fragments':len(frags),'junctions':len(jmap),'alternate_reads_with_CIGAR_N_inside_haplotype_window':allele_N}
expected_rows=list(csv.DictReader((W/'new-D8/junctions.tsv').open(),delimiter='\t'));mismatches=[]
for row in expected_rows:
 key=(row['gene'],int(row['intron_start0']),int(row['intron_end0']));v=observed_junctions.get(key)
 want={'all':int(row['query_name_fragments']),'alt':int(row['alternate_linked_pairs']),'direct':int(row['alternate_and_junction_same_read'])}
 if v!=want:mismatches.append({'key':key,'expected':want,'observed':v})
for s in S['summaries']:
 if s['allele_states']!=answer[s['gene']]['fragment_allele_states']:mismatches.append({'gene':s['gene'],'allele_count_mismatch':True})
report={'pysam_version':pysam.__version__,'independent_method':'Explicit CIGAR reference/query walk; no import of production audit; same stated filters','summary':answer,'junction_rows_compared':len(expected_rows),'independent_junction_rows':len(observed_junctions),'mismatches':mismatches,'primary_fragment_anomalies':odd_mates,'script_sha256':hashlib.sha256((W/'audit_target_splicing.py').read_bytes()).hexdigest()}
(W/'independent-new-D8-count-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

import pathlib,json,pysam,collections
from Bio import Align
O=pathlib.Path(__file__).parent;B=pathlib.Path.home()/'.local/share/codex/caris-analysis';N=json.loads((O/'nominations.json').read_text());fq=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'));rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
base=json.loads((B/'oct4-fusion-read-audit/patient/junction-read-audit.json').read_text());r=next(x for x in base['rows']if x['candidate_id']=='d97835c4fba3');query=next(x['sequence']for x in r['junction_patterns']if x['arm_bases']==25)
a=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-4)
for g in('insertion','deletion'):
 setattr(a,f'open_internal_{g}_score',-6);setattr(a,f'extend_internal_{g}_score',-1)
for s in('left','right'):
 setattr(a,f'open_{s}_deletion_score',0);setattr(a,f'extend_{s}_deletion_score',0);setattr(a,f'open_{s}_insertion_score',-2);setattr(a,f'extend_{s}_insertion_score',-2)
refs=collections.defaultdict(list);genomehits=set();txhits=set();nscored=0
# Every50nt query seed occurs among the already scanned pure-repeat control's phases. Reuse its nominated windows only; no new whole-genome scan.
for ori,qq in [('+',query),('-',rc(query))]:
 for st in sorted(set(list(range(0,len(qq)-18,19))+[len(qq)-19])):assert qq[st:st+19]in N['seed_patterns']
 for x in N['queries']:
  if x['query']['query_id']!='C_repeat_nonunique':continue
  for w in x['references']['genome']:
   ss=fq.fetch(w['reference_id'],w['start0'],w['end0']).upper();nscored+=1
   if a.score(ss,qq)==100:
    for j in range(len(ss)-49):
     if ss[j:j+50]==qq:genomehits.add((ori,w['reference_id'],w['start0']+j,w['start0']+j+50))
  for w in x['references']['transcript']:refs[w['reference_id']].append((ori,qq,w))
with pysam.FastxFile(str(B/'oct4-hla-allele-support/gencode.v37.transcripts.fa.gz'))as f:
 for t in f:
  tid=t.name.split('|')[0]
  for ori,qq,w in refs.get(tid,[]):
   ss=t.sequence[w['start0']:w['end0']].upper();nscored+=1
   if a.score(ss,qq)==100:
    for j in range(len(ss)-49):
     if ss[j:j+50]==qq:txhits.add((ori,tid,w['start0']+j,w['start0']+j+50))
assert len(genomehits)==14 and len(txhits)==10
report={'status':'pass','postrun_control_extension':'The161nt pure-repeat control had no perfect full-length placement. Added a known multiply represented50nt repeat to demonstrate exact-score ambiguity, without changing patient results or rerunning whole-reference seed nomination.','query':query,'expected_perfect_score':100,'genome_occurrences':len(genomehits),'transcript_occurrences':len(txhits),'genome_hits':sorted(genomehits),'transcript_hits':sorted(txhits),'nominated_windows_scored':nscored,'patient_results_unchanged':True}
(O/'short-repeat-control.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('genome_hits','transcript_hits')}))

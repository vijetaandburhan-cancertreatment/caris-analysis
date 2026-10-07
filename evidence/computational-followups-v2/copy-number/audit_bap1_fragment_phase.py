"""Independent bounded phase audit using pysam aligned-pair mapping.
Read-only resident DNA BAM; compact prior public UCSC JSON as reference.
"""
from pathlib import Path
from collections import defaultdict,Counter
import pysam,json,hashlib,statistics
ROOT=Path.cwd();O=ROOT/'work/oct4-followup/copy-number';R=json.loads((ROOT/'work/oct1-deep/splicing/reference/BAP1.hg38-reference.json').read_text());CH='chr3';SNP=52408337;DS=52408550;DE=52408561
ref=lambda a,b:R['dna'][a-R['start']:b-R['start']].upper()
assert ref(SNP,SNP+1)=='G';assert ref(DS-1,DE)=='CGCCGGGACCGG'
B=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/DNA_TN26-279853.bam');I=ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai'
with pysam.AlignmentFile(str(B),'rb',index_filename=str(I)) as bam:
 reads=[r for r in bam.fetch(CH,SNP-500,DE+500) if not r.flag&3852 and r.is_proper_pair and r.mapping_quality>=30 and r.query_qualities is not None]
results=[]
for flank in [3,10,20]:
 a=DS-flank;b=DE+flank;wt=ref(a,b);mut=ref(a,DS)+ref(DE,b);groups=defaultdict(list)
 for r in reads:
  mp={g:q for q,g in r.get_aligned_pairs(matches_only=False) if q is not None and g is not None};obs={};clean={};qS=mp.get(SNP)
  soft=any(op==4 for op,n in r.cigartuples)
  if qS is not None and r.query_qualities[qS]>=25 and r.query_sequence[qS] in 'ACGT':
   obs['SNP']=r.query_sequence[qS];clean['SNP']=(not soft and qS-r.query_alignment_start>=5 and r.query_alignment_end-qS-1>=5)
  if a in mp and b-1 in mp and not any(op==3 for op,n in r.cigartuples):
   q1,q2=mp[a],mp[b-1]+1
   if q2>q1 and min(r.query_qualities[q1:q2])>=25:
    h=r.query_sequence[q1:q2];obs['DEL']='del11' if h==mut else 'reference' if h==wt else 'other';clean['DEL']=(not soft and q1-r.query_alignment_start>=5 and r.query_alignment_end-q2>=5)
  groups[(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name)].append((r,obs,clean))
 counts=Counter();clean_counts=Counter();marg=Counter();details=[];issues=[];patternset=set()
 for key,entries in groups.items():
  sn={o['SNP'] for _,o,_ in entries if 'SNP' in o};de={o['DEL'] for _,o,_ in entries if 'DEL' in o}
  for locus,x in [('SNP',sn),('deletion',de)]:
   if len(x)==1:marg[locus+'_'+next(iter(x))+'_fragments']+=1
   if len(x)>1:marg[locus+'_discordant_fragments']+=1
  if not sn or not de:continue
  if len(sn)!=1 or len(de)!=1:counts['discordant']+=1;continue
  s,d=next(iter(sn)),next(iter(de));label=s+'__'+d;counts[label]+=1
  clean=all(any(o.get(l)==v and c.get(l) for _,o,c in entries) for l,v in [('SNP',s),('DEL',d)])
  if clean:clean_counts[label]+=1
  # Every alignment from the fragment, not just locus-supporting alignments.
  rr=[r for r,_,_ in entries];pairvalid=False;why=[]
  if len(rr)==2 and {r.is_read1 for r in rr}=={True,False}:
   x,y=rr;reciprocal=x.next_reference_id==y.reference_id and y.next_reference_id==x.reference_id and x.next_reference_start==y.reference_start and y.next_reference_start==x.reference_start
   tlen=x.template_length==-y.template_length and abs(x.template_length)==max(x.reference_end,y.reference_end)-min(x.reference_start,y.reference_start)
   opposing=x.is_reverse!=y.is_reverse
   pairvalid=reciprocal and tlen and opposing
   if not reciprocal:why.append('mate coordinates')
   if not tlen:why.append('TLEN')
   if not opposing:why.append('orientation')
  else:why.append('not exactly one primary R1/R2 pair')
  if not pairvalid:issues.append(dict(fragment=hashlib.sha256('|'.join(key).encode()).hexdigest(),issues=why))
  pat=tuple(sorted((r.reference_start,r.reference_end,r.cigarstring,r.is_read1,r.is_reverse,r.template_length) for r in rr));patternset.add(pat)
  details.append(dict(fragment_sha256=hashlib.sha256('|'.join(key).encode()).hexdigest(),phase=label,clean=clean,proper_reciprocal_pair_and_TLEN=pairvalid,same_read_phase=any('SNP' in o and 'DEL' in o for _,o,_ in entries),tlen=abs(rr[0].template_length),alignments=pat))
 results.append(dict(flank=flank,joint_counts=dict(counts),clean_counts=dict(clean_counts),marginal_counts=dict(marg),unique_joint_pair_patterns=len(patternset),joint_fragment_count=len(details),same_read_phase=sum(d['same_read_phase'] for d in details),pair_failures=issues,tlen_range=[min(d['tlen'] for d in details),max(d['tlen'] for d in details)],joint_details=details))
prior=json.loads((O/'BAP1-SNP-frameshift-fragment-phase.json').read_text())
for x,y in zip(results,prior['results']):
 assert x['flank']==y['flank_bases'];assert x['joint_counts']==y['joint_fragment_counts'];assert x['clean_counts']==y['joint_clean_fragment_counts'];assert not x['pair_failures']
 assert {d['fragment_sha256'] for d in x['joint_details']}=={d['fragment_identifier_sha256'] for d in y['joint_fragment_details']}
res=dict(method='Independent pysam get_aligned_pairs mapping, local sequence reconstruction, reciprocal mate and TLEN checks; reference from separately retained public BAP1 JSON. Same sample/BAM is not orthogonal clinical validation.',filters='exclude flags3852 (unmapped/mate-unmapped/secondary/QCfail/duplicate/supplementary); proper pair; MAPQ>=30; BQ>=25; exact haplotypes; discordant locus calls excluded; clean: no softclip and >=5 bases from aligned read ends',reference_source='work/oct1-deep/splicing/reference/BAP1.hg38-reference.json',reference_sha256=hashlib.sha256((ROOT/'work/oct1-deep/splicing/reference/BAP1.hg38-reference.json').read_bytes()).hexdigest(),passed=True,prior_counts_and_fragment_identifiers_agree=True,results=results)
(O/'BAP1-SNP-phase-independent-audit.json').write_text(json.dumps(res,indent=2)+'\n')
for x in results:print(x['flank'],x['joint_counts'],x['clean_counts'],'pair failures',len(x['pair_failures']),'patterns',x['unique_joint_pair_patterns'],'TLEN',x['tlen_range'])

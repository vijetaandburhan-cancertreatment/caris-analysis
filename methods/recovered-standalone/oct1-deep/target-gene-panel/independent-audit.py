"""Independent selected-transcript and bounded haplotype audit, no network."""
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import json,gzip,hashlib
import pysam
from Bio import SeqIO
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
RAW=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
A=json.loads((P/'target-annotation.json').read_text())['genes'];S=json.loads((P/'screen-results.json').read_text())
ALL={a['transcript']:a for a in json.loads((ROOT/'work/oct1-deep/coding-indels/selected-coding-transcripts.json').read_text())}
BASES='TCAG';AA='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';CODE=dict(zip((a+b+c for a in BASES for b in BASES for c in BASES),AA));COMP=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(COMP)[::-1]
def tr(s):return ''.join(CODE.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
proteins={}
with gzip.open(ROOT/'work/oct1-deep/neoantigen/gencode.v37.pc_translations.fa.gz','rt') as f:
 for r in SeqIO.parse(f,'fasta'):
  tx=r.id.split('|')[1]
  if any(a['selected_transcript']==tx for a in A.values()):proteins[tx]=str(r.seq)
anns=[];con=[]
for gene,a in A.items():
 tx=a['selected_transcript'];old=ALL[tx];w=json.loads((P/f'{gene}.hg38-reference.json').read_text());genomic=w['dna'].upper();pos=[p for s,e in sorted(a['selected_CDS']) for p in range(s,e)];nt=''.join(genomic[p-w['start']] for p in pos)
 if a['strand']=='-':nt=rc(nt);pos.reverse()
 idx={p:i for i,p in enumerate(pos)};prot=tr(nt)
 anns.append({'gene':gene,'transcript':tx,'strand':a['strand'],'MANE':'MANE_Select' in old['tags'],'tags':old['tags'],'selected_CDS_matches_source':sorted(a['selected_CDS'])==sorted(old['CDS']),'reference_protein_exact_match':prot==proteins[tx],'CDS_bases':len(nt)})
 for v in S['candidates']:
  if v['gene']!=gene:continue
  p=v['pos1']-1;ref=v['ref'];alt=v['alt'];r={'gene':gene,'chrom':v['chrom'],'pos1':v['pos1'],'ref':ref,'alt':alt,'reported_consequence':v['consequence'],'reported_protein':v['protein'],'selected_transcript':tx}
  while ref and alt and ref[0]==alt[0]:p+=1;ref=ref[1:];alt=alt[1:]
  while ref and alt and ref[-1]==alt[-1]:ref=ref[:-1];alt=alt[:-1]
  if (ref and not all(z in idx for z in range(p,p+len(ref)))) or (not ref and (p-1 not in idx or p not in idx)):
   r['independent_status']='noncoding_or_boundary';r['distance_to_nearest_CDS_base']=min(abs(v['pos1']-1-z) for z in pos);con.append(r);continue
  if ref:
   ids=sorted(idx[z] for z in range(p,p+len(ref)));assert ids==list(range(ids[0],ids[-1]+1));left,right=ids[0],ids[-1]+1
  else:
   ids=sorted([idx[p-1],idx[p]]);assert ids[1]-ids[0]==1;left=right=ids[1]
  rr=ref if a['strand']=='+' else rc(ref);aa=alt if a['strand']=='+' else rc(alt);assert nt[left:right]==rr
  mutant=nt[:left]+aa+nt[right:];mp=tr(mutant);pref=0
  while pref<min(len(prot),len(mp)) and prot[pref]==mp[pref]:pref+=1
  suffix=0
  while suffix<min(len(prot)-pref,len(mp)-pref) and prot[-suffix-1]==mp[-suffix-1]:suffix+=1
  endr=len(prot)-suffix if suffix else len(prot);endm=len(mp)-suffix if suffix else len(mp)
  first='=' if pref==len(prot)==len(mp) else prot[pref]+str(pref+1)+(mp[pref] if pref<len(mp) else '-')
  r.update(independent_status='single_edit_reference_translation',first_difference=first,first_difference_matches=first==v['protein'],coding_start1=left+1,coding_deleted=rr,coding_inserted=aa,net_base_change=len(aa)-len(rr),reference_protein_deleted=prot[pref:endr],mutant_protein_inserted=mp[pref:endm],left_context=prot[max(0,pref-10):pref],right_context=prot[endr:endr+10],reference_protein_length=len(prot),mutant_protein_length=len(mp));con.append(r)

def exact_read(read,start,end):
 rp=read.reference_start;qp=0;seq=[];qs=[];covers=set()
 for op,n in read.cigartuples:
  if op in (0,7,8):
   x,y=max(rp,start),min(rp+n,end)
   if x<y:
    i,j=qp+x-rp,qp+y-rp;seq.append(read.query_sequence[i:j]);qs.extend(read.query_qualities[i:j]);covers.update([x,y-1])
   rp+=n;qp+=n
  elif op==1:
   if start<rp<end:seq.append(read.query_sequence[qp:qp+n]);qs.extend(read.query_qualities[qp:qp+n])
   qp+=n
  elif op==2:rp+=n
  elif op==3:
   if rp<end and rp+n>start:return None,'splice_in_window'
   rp+=n
  elif op==4:qp+=n
 if start not in covers or end-1 not in covers:return None,'incomplete'
 if not qs or min(qs)<20:return None,'BQ'
 return ''.join(seq),'ok'

reads=[]
targets=[v for v in S['candidates'] if len(v['ref'])!=len(v['alt'])]
for kind in ['DNA','RNA']:
 with pysam.AlignmentFile(str(RAW/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai')) as bam:
  for v in targets:
   gene=v['gene'];w=json.loads((P/f'{gene}.hg38-reference.json').read_text());g=w['dna'].upper()
   for flank in [12,30]:
    pos=v['pos1']-1;start=pos-flank;end=pos+len(v['ref'])+flank;ref=g[start-w['start']:end-w['start']];assert ref[flank:flank+len(v['ref'])]==v['ref'];alt=ref[:flank]+v['alt']+ref[flank+len(v['ref']):]
    names=defaultdict(set);stats=Counter();evidence=defaultdict(list)
    for r in bam.fetch(v['chrom'],start,end):
     if r.flag&(4|256|512|1024|2048) or r.mapping_quality<20 or r.query_qualities is None:continue
     if r.has_tag('NH') and r.get_tag('NH')>1:continue
     seq,reason=exact_read(r,start,end)
     if seq is None:stats[reason]+=1;continue
     state='alt' if seq==alt else 'ref' if seq==ref else 'other';stats[state]+=1;name=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name);names[name].add(state)
     if state=='alt':evidence[name].append({'start0':r.reference_start,'end0':r.reference_end,'mate_start0':r.next_reference_start,'reverse':r.is_reverse,'cigar':r.cigarstring,'mapq':r.mapping_quality,'softclip_bases':sum(n for op,n in r.cigartuples if op==4),'nM':r.get_tag('nM') if r.has_tag('nM') else None})
    collapsed=Counter(next(iter(x)) if len(x)==1 else 'discordant' for x in names.values());altrows=[z for n,ls in evidence.items() if names[n]=={'alt'} for z in ls]
    reads.append({'gene':gene,'kind':kind,'chrom':v['chrom'],'pos1':v['pos1'],'ref':v['ref'],'alt':v['alt'],'flank':flank,'read_counts':dict(stats),'paired_names':dict(collapsed),'alt_forward_reads':sum(not r['reverse'] for r in altrows),'alt_reverse_reads':sum(r['reverse'] for r in altrows),'alt_unclipped_reads':sum(r['softclip_bases']==0 for r in altrows),'alt_evidence_first20':altrows[:20],'distinct_alt_endpoints':len({(r['start0'],r['end0'],r['reverse'],r['mate_start0']) for r in altrows})})
summary={'genes':len(A),'MANE':sum(r['MANE'] for r in anns),'candidates':len(S['candidates']),'consequences':dict(Counter(v['consequence'] for v in S['candidates'])),'CDS_plus4_bases':sum(g['CDS_plus4_bases'] for g in S['gene_summaries']),'bases_ge20':sum(g['bases_ge20'] for g in S['gene_summaries']),'all_reference_proteins_match':all(r['reference_protein_exact_match'] for r in anns),'all_CDS_protein_first_differences_match':all(r.get('first_difference_matches',True) for r in con)}
result={'created_utc':datetime.now(timezone.utc).isoformat(),'method':'Independent reference-coordinate transcript edit/translation and every selected reference protein compared with GENCODE37 FASTA. Bounded all-indel DNA/RNA exact12/30flank recount uses MAPQ20/BQ20, primary/nonduplicate/QCpass and NH1 where available, query names collapsed. No primary scan duplicate and no clinical classification.','summary':summary,'transcripts':anns,'consequences':con,'indel_haplotype_audit':reads}
(P/'independent-audit-data.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summary,indent=2))
for r in con:
 if r['gene']=='FGFR1' or r['independent_status']=='noncoding_or_boundary':print(json.dumps(r))
for r in reads:print(r['gene'],r['pos1'],r['ref'],r['alt'],r['kind'],r['flank'],r['paired_names'])

"""Bounded, local comparison of a discordant locus with its homologous pseudogene.
Only generic public reference windows are fetched; all reads stay local.
"""
from pathlib import Path
from collections import defaultdict,Counter
import json,gzip,re,urllib.request,datetime
import pysam
from Bio.Align import PairwiseAligner
P=Path(__file__).resolve().parent;ROOT=P.parents[2];RAW=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');genes={};exons=defaultdict(list)
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  a=line.rstrip().split('\t');m=re.search('gene_name "([^"]+)"',a[8])
  if not m or m[1] not in {'ZRSR2','ZRSR2P1'}:continue
  gene=m[1]
  if a[2]=='gene':genes[gene]={'chrom':a[0],'start':int(a[3])-1,'end':int(a[4]),'strand':a[6]}
  if a[2]=='exon' and ('transcript_id "ENST00000307771.8"' in a[8] or 'transcript_id "ENST00000512790.6"' in a[8]):exons[gene].append((int(a[3])-1,int(a[4])))
refs={};sequences={};position=15820242
for gene,g in genes.items():
 file=P/f'{gene}.public-reference.json';url=f'https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={g["chrom"]};start={g["start"]};end={g["end"]}'
 if not file.exists():
  with urllib.request.urlopen(url,timeout=45) as f:obj=json.load(f)
  obj['source_url']=url;obj['retrieved_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();file.write_text(json.dumps(obj)+'\n')
 w=json.loads(file.read_text());assert len(w['dna'])==g['end']-g['start'];refs[gene]=w['dna'].upper();sequences[gene]=''.join(refs[gene][s-g['start']:e-g['start']] for s,e in sorted(exons[gene]));assert g['strand']=='+'
positions=[p for s,e in sorted(exons['ZRSR2']) for p in range(s,e)];i=positions.index(position);xref=sequences['ZRSR2'];assert xref[i]=='C';xt=xref[:i]+'T'+xref[i+1:]
aligner=PairwiseAligner();aligner.mode='local';aligner.match_score=2;aligner.mismatch_score=-3;aligner.open_gap_score=-5;aligner.extend_gap_score=-1
rows=[];summary=[]
for kind in ['DNA','RNA']:
 with pysam.AlignmentFile(str(RAW/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai')) as bam:
  calls=defaultdict(set);details=[]
  for read in bam.fetch('chrX',position,position+1):
   if read.flag&(4|256|512|1024|2048) or read.mapping_quality<30 or read.query_qualities is None or (read.has_tag('NH') and read.get_tag('NH')!=1):continue
   rp=read.reference_start;qp=0;qsite=None
   for op,n in read.cigartuples:
    if op in (0,7,8):
     if rp<=position<rp+n:qsite=qp+position-rp
     qp+=n;rp+=n
    elif op in (1,4):qp+=n
    elif op in (2,3):rp+=n
   if qsite is None or read.query_qualities[qsite]<30 or min(qsite-read.query_alignment_start,read.query_alignment_end-qsite-1)<5:continue
   base=read.query_sequence[qsite];name=(read.get_tag('RG') if read.has_tag('RG') else '')+'|'+read.query_name;calls[name].add(base)
   if kind=='RNA':
    query=read.query_alignment_sequence;scores={'ZRSR2_reference_C':aligner.score(xref,query),'ZRSR2_observed_DNA_T':aligner.score(xt,query),'ZRSR2P1_reference':aligner.score(sequences['ZRSR2P1'],query)}
    detail={'name':name,'base':base,'mapq':read.mapping_quality,'NH':read.get_tag('NH') if read.has_tag('NH') else None,'start0':read.reference_start,'end0':read.reference_end,'cigar':read.cigarstring,'reverse':read.is_reverse,'aligned_query_length':len(query),'alignment_scores':scores,'difference_P1_minus_best_X':scores['ZRSR2P1_reference']-max(scores['ZRSR2_reference_C'],scores['ZRSR2_observed_DNA_T']),'query_sequence_local_only':query};details.append(detail)
  count=Counter(next(iter(v)) if len(v)==1 else 'discordant' for v in calls.values());summary.append({'kind':kind,'counts':dict(count)})
  if kind=='RNA':rows=details
scoreclasses=defaultdict(Counter)
for r in rows:
 scoreclasses[r['base']]['P1_better' if r['difference_P1_minus_best_X']>0 else 'tie' if r['difference_P1_minus_best_X']==0 else 'X_better']+=1
out={'scope':'Independent CIGAR-based recount of single flagged locus. Local alignment score comparison of the aligned portion of each qualifying RNA read to GENCODE37 ZRSR2 full transcript (reference C or observed DNA T at the flagged site) versus reference ZRSR2P1 transcript. Generic public references only; no external patient reads. Scores use +2 match, -3 mismatch, -5 gap open, -1 extension; this exploratory comparison is not genome-wide remapping or calibrated likelihood.','public_reference_genes':genes,'independent_counts':summary,'RNA_read_score_classes':{k:dict(v) for k,v in scoreclasses.items()},'RNA_reads':rows,'limitations':['Comparison to two reference transcripts does not prove a read source or explain all ambiguity; actual haplotypes and other homologous loci may differ.','Read-level score counts are not independent fragments or a contamination estimate.','One discordant locus does not establish sample mismatch, contamination, RNA editing or pathology.']}
(P/'audit-zrsr2.json').write_text(json.dumps(out,indent=2)+'\n');print(summary);print(out['RNA_read_score_classes'])

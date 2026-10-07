"""Read-only bounded splice inventory from existing STAR BAM; research, not clinical calling."""
from pathlib import Path
import collections,csv,gzip,hashlib,json,re,resource,time
import pysam
OUT=Path(__file__).resolve().parent; ROOT=OUT.parents[2]
GTF=ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz'
BAM=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/RNA_TN26-279853.bam')
INDEX=ROOT/'work/oct1-analysis/RNA_TN26-279853.bam.bai'
PANEL=set('BAP1 RASA1 APC LATS1 LATS2 NF2 SMARCA4 SMARCB1 SMARCA2 TP53 CDKN2A CDKN2B MTAP SETD2 PBRM1 RB1 PTEN ARID1A ARID1B STAG2 KDM6A TRAF7 SUFU VHL B2M JAK1 JAK2 MET ALK ROS1 RET NTRK1 NTRK2 NTRK3 FGFR1 FGFR2 FGFR3 BRAF RAF1 EGFR ERBB2 ERBB3 ERBB4 KIT PDGFRA PDGFRB NRG1 AXL MAP2K1 MAP2K2 AKT1 AKT2 AKT3 PIK3CA KRAS NRAS HRAS'.split())
# Canonical chromosomes match the existing GTF; no inferred strands from read orientation.
start_time=time.time(); stats=collections.Counter(); all_tx={}; targets={}; target_tx={}
print('Parsing GENCODE37 annotation',flush=True)
with gzip.open(GTF,'rt') as fh:
 for ln in fh:
  stats['gtf_lines']+=1
  if ln.startswith('#'):continue
  a=ln.rstrip().split('\t')
  if a[2] not in ('gene','exon','transcript','CDS'):continue
  gene_match=re.search(r'gene_name "([^"]+)"',a[8]); gene=gene_match.group(1) if gene_match else None
  if a[2]!='exon' and gene not in PANEL:continue
  attrs=collections.defaultdict(list)
  for k,v in re.findall(r'(\w+) "([^"]+)"',a[8]): attrs[k].append(v)
  chrom=a[0]; s=int(a[3])-1; e=int(a[4]); strand=a[6]
  if a[2]=='gene' and gene in PANEL:
   targets[gene]={'gene':gene,'gene_id':attrs['gene_id'][0],'chrom':chrom,'start0':s,'end0':e,'strand':strand}
  if 'transcript_id' not in attrs:continue
  tid=attrs['transcript_id'][0]
  if a[2]=='exon':
   tx=all_tx.setdefault(tid,[chrom,[]]); tx[1].append((s,e))
  if gene in PANEL:
   tx=target_tx.setdefault(tid,{'transcript_id':tid,'gene':gene,'chrom':chrom,'strand':strand,'exons':[],'CDS':[],'tags':[]})
   tx['tags']=sorted(set(tx['tags']+attrs['tag']))
   if a[2]=='exon':tx['exons'].append([s,e])
   if a[2]=='CDS':tx['CDS'].append([s,e])
annotated=set()
for tid,(chrom,exons) in all_tx.items():
 exons.sort()
 for (s1,e1),(s2,e2) in zip(exons,exons[1:]):
  if s2>e1: annotated.add((chrom,e1,s2))
stats['annotated_unique_junctions']=len(annotated); stats['annotation_transcripts']=len(all_tx)
del all_tx
known_target=collections.defaultdict(list); target_sites=collections.defaultdict(lambda:collections.defaultdict(set)); skipped=collections.defaultdict(list)
for tid,tx in target_tx.items():
 exons=sorted(tx['exons']); tx['exons']=exons; tx['CDS'].sort()
 for i,((s1,e1),(s2,e2)) in enumerate(zip(exons,exons[1:])):
  key=(tx['chrom'],e1,s2);known_target[key].append(tid)
  target_sites[tx['gene']]['left'].add(e1);target_sites[tx['gene']]['right'].add(s2)
 # Every exon-skipping junction in an explicitly annotated transcript, capped to five skipped exons.
 for i in range(len(exons)):
  for j in range(i+2,min(len(exons),i+7)):
   omitted=exons[i+1:j]; coding=sum(max(0,min(ee,ce)-max(ss,cs)) for ss,ee in omitted for cs,ce in tx['CDS'])
   skipped[(tx['chrom'],exons[i][1],exons[j][0])].append({'transcript':tid,'gene':tx['gene'],'skipped_genomic_exons_1based':list(range(i+2,j+1)),'skipped_transcript_exons_1based':list(range(i+2,j+1)) if tx['strand']=='+' else list(range(len(exons)-j+1,len(exons)-i)), 'skipped_CDS_nt':coding,'frame_mod3':coding%3})

def priority(tx):
 tags=tx['tags']; appris=[int(t.rsplit('_',1)[1]) for t in tags if t.startswith('appris_principal_')]
 return ('MANE_Select' in tags,-min(appris) if appris else -99,'basic' in tags,sum(e-s for s,e in tx['CDS']))
for gene,g in targets.items():
 ts=[t for t in target_tx.values() if t['gene']==gene and t['CDS']]
 g['selected_transcript']=max(ts,key=priority)['transcript_id'] if ts else None
by_chrom=collections.defaultdict(list)
for gene,g in targets.items(): by_chrom[g['chrom']].append(g)
annotation={'GTF':str(GTF),'GTF_sha256':hashlib.sha256(GTF.read_bytes()).hexdigest(),'coordinate_system':'0-based half-open intron [start,end), connecting exonic base start-1 to end','genes':targets,'transcripts':target_tx,'global_annotated_junctions':len(annotated),'panel_genes':sorted(PANEL)}
(OUT/'target-annotation.json').write_text(json.dumps(annotation,indent=2)+'\n')
print('Annotation ready:',len(annotated),'junctions;',len(targets),'target genes; elapsed',round(time.time()-start_time,1),flush=True)
# Genome-wide counts are alignment/read counts. Target counts additionally collapse RG+query name.
counts={}; detailed={}
MASK=4|256|512|1024|2048; MATCH={0,7,8}; CONSUME_R={0,2,3,7,8}; CONSUME_Q={0,1,4,7,8}
with pysam.AlignmentFile(str(BAM),'rb',index_filename=str(INDEX),threads=1) as bam:
 for read in bam.fetch(until_eof=True):
  stats['BAM_records']+=1
  if stats['BAM_records']%5000000==0:
   print('Read',stats['BAM_records'],'junctions',len(counts),'target junctions',len(detailed),'seconds',round(time.time()-start_time,1),'RSS_MB',round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,1),flush=True)
  if read.flag&MASK:continue
  stats['eligible_primary_flags']+=1
  if read.mapping_quality<20:continue
  if read.has_tag('NH') and read.get_tag('NH')!=1:continue
  stats['MQ20_unique_primary']+=1
  cigar=read.cigartuples
  if not cigar or not any(op==3 for op,n in cigar):continue
  stats['spliced_primary_reads_before_anchor']+=1
  chrom=bam.references[read.reference_id]; rpos=read.reference_start;qpos=0
  qual=read.query_qualities; xs=read.get_tag('XS') if read.has_tag('XS') else '.'
  nm=read.get_tag('nM') if read.has_tag('nM') else (read.get_tag('NM') if read.has_tag('NM') else None)
  for ci,(op,n) in enumerate(cigar):
   if op==3:
    stats['CIGAR_N_events']+=1
    if not (20<=n<=500000):stats['intron_length_failed']+=1
    elif ci==0 or ci==len(cigar)-1 or cigar[ci-1][0] not in MATCH or cigar[ci+1][0] not in MATCH:stats['adjacent_match_failed']+=1
    else:
     left,right=cigar[ci-1][1],cigar[ci+1][1]
     if min(left,right)<12:stats['anchor12_failed']+=1
     elif qual is None or min(qual[qpos-12:qpos])<20 or min(qual[qpos:qpos+12])<20:stats['BQ20_anchor_failed']+=1
     else:
      key=(chrom,rpos,rpos+n); stats['qualified_junction_read_events']+=1
      strict20=min(left,right)>=20 and min(qual[qpos-20:qpos])>=20 and min(qual[qpos:qpos+20])>=20
      # reads, 20bp reads, forward, reverse, zero mismatch, plus/minus/unknown XS, largest joint anchor
      vals=counts.setdefault(key,[0]*9);vals[0]+=1;vals[1]+=strict20;vals[3 if read.is_reverse else 2]+=1;vals[4]+=nm==0;vals[5 if xs=='+' else 6 if xs=='-' else 7]+=1;vals[8]=max(vals[8],min(left,right))
      genes=[g['gene'] for g in by_chrom[chrom] if g['start0']<=rpos and rpos+n<=g['end0']]
      if genes:
       d=detailed.setdefault(key,{'genes':genes,'fragments':set(),'strict20_fragments':set(),'signatures':set(),'zero_mismatch_fragments':set(),'examples':[],'reverse_fragments':set(),'forward_fragments':set()})
       fragment=(read.get_tag('RG') if read.has_tag('RG') else '')+'|'+read.query_name
       d['fragments'].add(fragment)
       if strict20:d['strict20_fragments'].add(fragment)
       if nm==0:d['zero_mismatch_fragments'].add(fragment)
       d['signatures'].add((read.reference_start,read.reference_end,read.is_reverse))
       d['reverse_fragments' if read.is_reverse else 'forward_fragments'].add(fragment)
       if len(d['examples'])<12:
        d['examples'].append({'query_name':read.query_name,'RG':read.get_tag('RG') if read.has_tag('RG') else '', 'start0':read.reference_start,'end0':read.reference_end,'CIGAR':read.cigarstring,'flag':read.flag,'MAPQ':read.mapping_quality,'NH':read.get_tag('NH') if read.has_tag('NH') else None,'nM_or_NM':nm,'XS':xs,'left_anchor':left,'right_anchor':right,'left_12bp':read.query_sequence[qpos-12:qpos],'right_12bp':read.query_sequence[qpos:qpos+12],'left_12_BQ':list(qual[qpos-12:qpos]),'right_12_BQ':list(qual[qpos:qpos+12])})
   if op in CONSUME_R:rpos+=n
   if op in CONSUME_Q:qpos+=n

headers=['chrom','intron_start0','intron_end0','annotated_GENCODE37','read_events','anchor20_read_events','forward_read_events','reverse_read_events','zero_mismatch_read_events','XS_plus','XS_minus','XS_missing','max_min_anchor']
with gzip.open(OUT/'all-junctions.tsv.gz','wt') as fh:
 w=csv.writer(fh,delimiter='\t');w.writerow(headers)
 for key,vals in sorted(counts.items()):w.writerow([*key,int(key in annotated),*vals])
rows=[]; examples={}
for key,d in detailed.items():
 vals=counts[key];chrom,s,e=key
 for gene in d['genes']:
  sites=target_sites[gene]; ann=key in annotated; gene_transcripts=[tid for tid in known_target[key] if target_tx[tid]['gene']==gene]
  skips=[z for z in skipped[key] if z['gene']==gene]
  typ='annotated_target_gene' if gene_transcripts else 'annotated_other_gene' if ann else 'unannotated_exon_skip' if skips else 'unannotated_pair_known_sites' if s in sites['left'] and e in sites['right'] else 'unannotated_one_known_site' if s in sites['left'] or e in sites['right'] else 'unannotated_both_sites'
  candidate=not ann and len(d['fragments'])>=5 and len(d['signatures'])>=3 and len(d['strict20_fragments'])>=3
  row={'gene':gene,'chrom':chrom,'intron_start0':s,'intron_end0':e,'intron_length':e-s,'classification':typ,'annotated_GENCODE37':int(ann),'fragments':len(d['fragments']),'anchor20_fragments':len(d['strict20_fragments']),'distinct_alignment_signatures':len(d['signatures']),'zero_mismatch_fragments':len(d['zero_mismatch_fragments']),'forward_fragments':len(d['forward_fragments']),'reverse_fragments':len(d['reverse_fragments']),'reads':vals[0],'XS_plus':vals[5],'XS_minus':vals[6],'XS_missing':vals[7],'target_strand':targets[gene]['strand'],'max_min_anchor':vals[8],'selected_transcript_annotated':int(targets[gene]['selected_transcript'] in gene_transcripts),'annotated_transcripts':';'.join(gene_transcripts),'exon_skip_descriptions':json.dumps(skips,separators=(',',':')),'passes_review_priority':int(candidate)}
  rows.append(row)
  if candidate or gene=='MET':examples[f'{gene}:{chrom}:{s}-{e}']={'row':row,'read_examples':d['examples']}
# Local donor/acceptor fractions are descriptive, not transcript PSI.
for row in rows:
 same_left=[x for x in rows if x['gene']==row['gene'] and x['intron_start0']==row['intron_start0']]
 same_right=[x for x in rows if x['gene']==row['gene'] and x['intron_end0']==row['intron_end0']]
 row['fraction_same_left_boundary']=round(row['fragments']/sum(x['fragments'] for x in same_left),6)
 row['fraction_same_right_boundary']=round(row['fragments']/sum(x['fragments'] for x in same_right),6)
with (OUT/'target-junctions.tsv').open('w') as fh:
 w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(sorted(rows,key=lambda x:(x['gene'],x['intron_start0'],x['intron_end0'])))
(OUT/'candidate-read-evidence.json').write_text(json.dumps(examples,indent=2)+'\n')
summary={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'elapsed_seconds':time.time()-start_time,'peak_RSS_MB_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,'pysam_version':pysam.__version__,'BAM':str(BAM),'BAM_bytes':BAM.stat().st_size,'GTF_sha256':annotation['GTF_sha256'],'stats':dict(stats),'all_junctions':len(counts),'all_annotated_junctions':sum(k in annotated for k in counts),'target_junction_gene_rows':len(rows),'target_review_candidates':sum(r['passes_review_priority'] for r in rows),'method':'Primary mapped; excludes secondary/supplementary/QCfail/duplicate flags; MAPQ>=20; NH=1 if present; intron20..500000; immediately adjacent matched CIGAR anchors>=12bp each with BQ>=20. Strict-anchor count additionally requires20bp each BQ>=20. Target counts collapse RG+queryname per junction; signatures are alignment start/end/strand, not UMIs. Global inventory counts read events. GENCODE annotation membership is coordinate-only; XS intron-motif strand counts retained, not inferred library strandedness. Candidate prioritization: unannotated globally, >=5 fragment names, >=3 start/end/strand signatures, >=3 strict20 fragment names. No matched normal or calibrated PSI/clinical validation.', 'limits':['BAM has no RNA duplicate marks: fragment names do not prove independent PCR molecules.','No reference mismatch checks applied to anchors in inventory; nM=0 counts retained for stricter prioritization.','Annotations do not capture all normal isoforms; absent annotation is not tumor specificity.','Existing STAR alignment max intron200000 means this inventory is not sensitive to all longer events.','No chimeric-junction output or supplementary alignments; this is not a full fusion screen.','Short reads cannot reconstruct full isoforms or establish peptide expression.','Tumor/normal/brain/immune-cell admixture and technical artifacts may produce apparent events.'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'inventory-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)

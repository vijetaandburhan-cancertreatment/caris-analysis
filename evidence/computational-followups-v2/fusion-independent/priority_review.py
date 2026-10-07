#!/usr/bin/env python3
"""Independent bounded four-candidate exact-sequence and DNA-locus audit."""
import json,gzip,collections,hashlib,pathlib,re,itertools
import pysam
HOME=pathlib.Path.home();BASE=HOME/'.local/share/codex/caris-analysis'
OUT=pathlib.Path(__file__).parent;DUR=BASE/'oct5-fusion-priority-review';DUR.mkdir(exist_ok=True)
AUD=BASE/'oct4-fusion-read-audit';data=json.loads((AUD/'patient/junction-read-audit.json').read_text());refs=json.loads((AUD/'patient-reference/reference-ambiguity.json').read_text())
ids=['078bc453229c','e2a8e3625460','d97835c4fba3','cb8931f687fc']
rows={r['candidate_id']:r for r in data['rows'] if r['candidate_id'] in ids}
fa=pysam.FastaFile(str(BASE/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'))
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def h(s):return hashlib.sha256(s.encode()).hexdigest()
def norm(n):return n[:-2] if n.endswith(('/1','/2')) else n
result={'status':'complete','source':'independent exact-sequence implementation; selected RNA FASTQ audit plus indexed local DNA BAM fetch only','candidates':{},'sources':{'arriba_output_format':'https://github.com/suhrig/arriba/wiki/05-Output-files','raw_audit':str(AUD/'patient/junction-read-audit.json'),'reference_audit':str(AUD/'patient-reference/reference-ambiguity.json')}}
markers={}
for cid,row in rows.items():
 seq=next(p['sequence'] for p in row['junction_patterns'] if p['arm_bases']==25)
 assert len(seq)==50
 markers[cid]=seq
 keep=['gene1','gene2','type','confidence','breakpoint1','breakpoint2','site1','site2','reading_frame','strand1(gene/fusion)','strand2(gene/fusion)','transcript_id1','transcript_id2','peptide_sequence','fusion_transcript','distinct_named_reads']
 d={k:row.get(k) for k in keep};d['junction_25_plus_25']=seq[:25]+'|'+seq[25:];d['reference_matches']=[]
 for i,p in enumerate(refs['patterns']):
  roles=[x for x in p['roles'] if x['candidate']==cid and x['arm']==25]
  if roles:
   idx=str(i);d['reference_matches'].append({'pattern_index':i,'sequence':p['sequence'],'roles':roles,'genome_occurrences':refs['genome']['counts'].get(idx,0),'transcript_occurrences':refs['transcripts']['occurrences'].get(idx,0),'transcript_genes':refs['transcripts']['matching_genes'].get(idx,{}),'genome_first_loci':refs['genome']['first_loci'].get(idx,[])[:6]})
 if row['gene1']==row['gene2']:
  chrom,sp1=row['breakpoint1'].split(':');_,sp2=row['breakpoint2'].split(':');p1,p2=int(sp1),int(sp2)
  refseq=fa.fetch(chrom,p1-25,p1).upper()+fa.fetch(chrom,p2-1,p2+24).upper()
  d['reference_breakend_marker']=refseq[:25]+'|'+refseq[25:]
  d['reference_mismatch']=[{'marker_position1':i+1,'reference':a,'observed':b,'genomic_position1':p1-25+i+1 if i<25 else p2+i-25} for i,(a,b) in enumerate(zip(refseq,seq)) if a!=b]
  d['coordinate_duplication_span_bp']=p1-p2+1
  if p1-p2<100:d['reference_duplicated_segment']=fa.fetch(chrom,p2-1,p1).upper()
 result['candidates'][cid]=d
# Independent scanner, no automaton or shared parser implementation.
stats={ci:{q:{'pairs':set(),'families':set(),'R1':0,'R2':0,'both_mates':0}for q in (20,30)}for ci in ids};n=0
with pysam.FastxFile(str(AUD/'patient/audit-selected.R1.fastq.gz')) as a,pysam.FastxFile(str(AUD/'patient/audit-selected.R2.fastq.gz')) as b:
 for r1,r2 in itertools.zip_longest(a,b):
  assert r1 and r2 and norm(r1.name)==norm(r2.name);n+=1
  ss=[r1.sequence,r2.sequence];qq=[r1.quality,r2.quality];fam=min(ss[0]+'|'+ss[1],rc(ss[1])+'|'+rc(ss[0]))
  for ci,m in markers.items():
   per=[[],[]]
   for j in range(2):
    for word in (m,rc(m)):
     k=ss[j].find(word)
     while k>=0:
      per[j].append(min(ord(x)-33 for x in qq[j][k:k+len(word)]));k=ss[j].find(word,k+1)
   for q in (20,30):
    hit=[any(v>=q for v in p)for p in per]
    if any(hit):
     z=stats[ci][q];z['pairs'].add(norm(r1.name));z['families'].add(h(fam));z['R1']+=hit[0];z['R2']+=hit[1];z['both_mates']+=all(hit)
for ci in ids:
 result['candidates'][ci]['independent_raw_50mer']={q:{k:len(v) if isinstance(v,set) else v for k,v in x.items()}for q,x in stats[ci].items()}
 for q in (20,30):
  expected=data['statistics'][ci]['raw_exact']['25'][f'Q{q}'];x=result['candidates'][ci]['independent_raw_50mer'][q]
  assert x['pairs']==expected['pairs'] and x['families']==expected['orientation_normalized_sequence_families'] and all(x[k]==expected[k]for k in ('R1','R2','both_mates'))
result['selected_RNA_pairs_checked']=n;result['RNA_recount_all_8_quality_counts_pass']=True
# DNA: fetch narrow union covering each relevant breakend. No whole-DNA scan.
bam=pysam.AlignmentFile(str(BASE/'TN26-279853/DNA_TN26-279853.bam'),'rb',index_filename=str(pathlib.Path.cwd()/'work/oct1-analysis/DNA_TN26-279853.bam.bai'))
for ci in ids[:2]:
 row=rows[ci];chrom,sp1=row['breakpoint1'].split(':');_,sp2=row['breakpoint2'].split(':');p1,p2=int(sp1),int(sp2)
 intervals=[(p2-501,p1+500)]if p1-p2<1000 else[(p1-501,p1+500),(p2-501,p2+500)]
 seen=set();records=[]
 for s,e in intervals:
  for a in bam.fetch(chrom,s,e):
   key=(a.query_name,a.flag,a.reference_start,a.cigarstring)
   if key not in seen:seen.add(key);records.append(a)
 strict=[a for a in records if not(a.flag&3844) and a.mapping_quality>=20 and a.query_sequence and a.query_qualities]
 dna={'fetch_intervals0':intervals,'alignments_all':len(records),'alignments_primary_MAPQ20_nonduplicate_QCpass':len(strict),'unique_fragments_strict':len({(a.query_name,a.get_tag('RG')if a.has_tag('RG')else'')for a in strict}),'exact_markers':{},'breakend_aligned_coverage_Q20':{},'nearby_CIGAR_insertions':[]}
 for label,word in [('observed_RNA_junction',markers[ci]),('reference_breakend_junction',result['candidates'][ci]['reference_breakend_marker'].replace('|',''))]:
  for minq in (20,30):
   supp={};examples=[]
   for a in strict:
    sequence=a.query_sequence;qualities=a.query_qualities
    for w in (word,rc(word)):
     start=sequence.find(w)
     while start>=0:
      if min(qualities[start:start+len(w)])>=minq:
       key=(a.query_name,a.get_tag('RG')if a.has_tag('RG')else'');supp[key]=a
       if len(examples)<8:examples.append({'name_sha256':h(a.query_name),'flag':a.flag,'pos1':a.reference_start+1,'cigar':a.cigarstring,'MAPQ':a.mapping_quality,'marker_start0':start,'SA':a.get_tag('SA')if a.has_tag('SA')else None})
      start=sequence.find(w,start+1)
   dna['exact_markers'][label+'_Q'+str(minq)]={'fragments':len(supp),'examples':examples}
 for p in (p1,p2):
  names=set();bases=collections.defaultdict(set)
  for a in strict:
   for qp,rp in a.get_aligned_pairs(matches_only=True):
    if rp==p-1 and a.query_qualities[qp]>=20:
     key=(a.query_name,a.get_tag('RG')if a.has_tag('RG')else'');names.add(key);bases[a.query_sequence[qp]].add(key)
  dna['breakend_aligned_coverage_Q20'][str(p)]={'fragments':len(names),'bases':{k:len(v)for k,v in bases.items()}}
 # Enumerate all >=3bp insertions within100bp of either breakend, then exact10bp flanks against genome.
 insertions=collections.defaultdict(lambda:{'reads':set(),'strict_10bp_flank_fragments':set(),'cigars':collections.Counter(),'strand':collections.Counter()})
 for a in strict:
  qp=0;rp=a.reference_start
  for op,l in a.cigartuples:
   if op==1 and l>=3 and min(abs(rp-p1),abs(rp-(p2-1)))<=100:
    ins=a.query_sequence[qp:qp+l];z=insertions[(rp,ins)];name=(a.query_name,a.get_tag('RG')if a.has_tag('RG')else'');z['reads'].add(name);z['cigars'][a.cigarstring]+=1;z['strand']['reverse'if a.is_reverse else'forward']+=1
    if qp>=10 and qp+l+10<=a.query_length:
     hap=a.query_sequence[qp-10:qp+l+10]
     if hap==fa.fetch(chrom,rp-10,rp).upper()+ins+fa.fetch(chrom,rp,rp+10).upper() and min(a.query_qualities[qp-10:qp+l+10])>=20:z['strict_10bp_flank_fragments'].add(name)
   if op in (0,7,8):qp+=l;rp+=l
   elif op in (1,4):qp+=l
   elif op in (2,3):rp+=l
 for (pos,ins),z in sorted(insertions.items(),key=lambda kv:-len(kv[1]['reads'])):
  dna['nearby_CIGAR_insertions'].append({'insertion_after_pos1':pos,'inserted_sequence':ins,'length':len(ins),'fragments':len(z['reads']),'strict10bpBQ20_flank_fragments':len(z['strict_10bp_flank_fragments']),'strand_read_alignments':dict(z['strand']),'top_cigars':z['cigars'].most_common(4)})
 result['candidates'][ci]['bounded_DNA']=dna
fa.close();bam.close()
result['limits']=['Exact markers are based on caller consensus; matching raw sequence is corroboration not a validated genomic rearrangement.','Reference exact-match counts do not cover near-matching paralogs or unrepresented normal haplotypes.','Read names and sequence families are not UMI-defined molecules. No matched-normal sample was used; no somatic or tumor-cell origin inference.','DNA check uses indexed reads aligned near the named breakends, not all unmapped DNA reads or a genome-wide structural-variant call; lack of local support is not exclusion.','Transcript frame refers to the supplied GENCODE37 isoform/model; no protein production or actionability is established.']
for folder in (OUT,DUR):(folder/'priority-review.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'complete','selected_RNA_pairs':n,'candidates':{ci:{'gene':x['gene1'],'RNA_recount':x['independent_raw_50mer'],'DNA':x.get('bounded_DNA')}for ci,x in result['candidates'].items()}}))

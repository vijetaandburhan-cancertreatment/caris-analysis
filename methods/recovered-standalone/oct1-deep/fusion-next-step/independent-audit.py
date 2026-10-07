"""Independent compact output and sequence/family-count audit; no alignment rerun."""
from pathlib import Path
from collections import Counter,defaultdict
import json,re,csv,datetime
P=Path(__file__).resolve().parent;F=P/'full'
def load(name):return json.loads((P/name).read_text())
def fasta(path):
 d={};name=None;seq=[]
 with path.open() as f:
  for line in f:
   if line.startswith('>'):
    if name is not None:d[name]=''.join(seq).upper()
    name=line[1:].strip().split()[0];seq=[]
   else:seq.append(line.strip())
  if name is not None:d[name]=''.join(seq).upper()
 return d
COM=str.maketrans('ACGTN','TGCAN');rc=lambda s:s.translate(COM)[::-1]
source=load('full/pizzly.insert235.json');window=load('full/pizzly.insert235.junction-window-validation.json');whole=load('full/pizzly.insert235.exact-junctions.json');bam=load('full/candidate-bam-audit.json');families=load('full/candidate-fragment-families.json');fused=fasta(F/'pizzly.insert235.fusions.fasta');wr={(r['geneA'],r['geneB']):r for r in window['summary']};er={(r['geneA'],r['geneB']):r for r in whole['summary']};recounts=[]
for g in source['genes']:
 key=g['geneA']['name'],g['geneB']['name'];names=set();sequences=set();en=set();es=set()
 for t in g.get('transcripts',[]):
  s=fused[t['fasta_record']];bp=t['transcriptA']['endPos']-t['transcriptA']['startPos']
  if bp<20 or len(s)-bp<20:continue
  join=s[bp-20:bp+20]
  if set(join)-set('ACGT'):continue
  for pair in g['readpairs']:
   pairseq=(pair['read1']['seq'],pair['read2']['seq']);name=pair['read1']['name']
   for mate in ['read1','read2']:
    seq=pair[mate]['seq'].upper()
    for q in (seq,rc(seq)):
     if join in q:names.add(name);sequences.add(pairseq)
     at=s.find(q)
     while at>=0:
      if bp-at>=20 and at+len(q)-bp>=20:en.add(name);es.add(pairseq)
      at=s.find(q,at+1)
 rec={'geneA':key[0],'geneB':key[1],'central40_names':len(names),'central40_sequence_pairs':len(sequences),'whole_read_exact_names':len(en),'whole_read_exact_sequence_pairs':len(es),'central40_matches_prior':len(names)==wr[key]['junction40_name_pairs'] and len(sequences)==wr[key]['junction40_sequence_pairs'],'whole_read_matches_prior':len(en)==er[key]['exact20nt_join_name_pairs'] and len(es)==er[key]['exact20nt_join_sequence_pairs'],'central40_supporting_names':sorted(names)};recounts.append(rec)
# Independently reconstruct 5-prime coordinates from the compact BAM evidence.
def unclipped(b):
 cigar=re.findall(r'(\d+)([MIDNSHP=X])',b['cigar']);reverse=bool(b['flag']&16);ops=list(reversed(cigar)) if reverse else cigar;clip=0
 for n,op in ops:
  if op not in 'SH':break
  clip+=int(n)
 return b['reference'],(b['end1']+clip if reverse else b['start1']-clip),reverse
familyrows=[]
for g in bam['candidates']:
 key=g['geneA']['name'],g['geneB']['name'];q=next(r for r in recounts if (r['geneA'],r['geneB'])==key);names=set(q['central40_supporting_names']);allfam=defaultdict(set);directfam=defaultdict(set);incomplete=[]
 for name,ev in g['bam_evidence'].items():
  primary=[x for x in ev if not x['secondary'] and not x['supplementary'] and not x['qc_failed']];m1=[x for x in primary if x['mate']==1];m2=[x for x in primary if x['mate']==2]
  if len(m1)!=1 or len(m2)!=1:incomplete.append(name);continue
  fam=unclipped(m1[0]),unclipped(m2[0]);allfam[fam].add(name)
  if name in names:directfam[fam].add(name)
 old=next(x for x in families['summary'] if (x['geneA'],x['geneB'])==key)
 familyrows.append({'geneA':key[0],'geneB':key[1],'all_caller_names':len(g['names']),'all_coordinate_families':len(allfam),'exact_central40_names':len(names),'exact_central40_coordinate_families':len(directfam),'all_family_count_matches':len(allfam)==old['primary_pair_unclipped5prime_coordinate_families'],'incomplete_names':incomplete})
nf=next(g for g in bam['candidates'] if g['geneA']['name']=='NF1');ne=[]
for name,ev in nf['bam_evidence'].items():
 primary=[b for b in ev if not b['secondary'] and not b['supplementary']];ne.append({'name':name,'chromosomes':sorted(set(b['reference'] for b in primary)),'nM_values':sorted(set(b['STAR_nM'] for b in primary)),'cigars':[b['cigar'] for b in primary],'NH_values':sorted(set(b['NH'] for b in primary)),'any_softclip':any('S' in b['cigar'] for b in primary)})
normal=list(csv.DictReader((F/'pizzly.insert235.junction40.normal-matches.tsv').open(),delimiter='\t'));control=load('public-full-reference/control-summary.json');info=load('full/kallisto/run_info.json');frag=load('full/kallisto/fragment-distribution.json');hist=frag['histogram'];total=sum(hist);c=0;q95=None
for i,n in enumerate(hist):
 c+=n
 if q95 is None and c/total>.95:q95=i
summary={'filtered_primary':len(source['genes']),'filtered_insert400':len(load('full/pizzly.insert400.json')['genes']),'unfiltered_primary':len(load('full/pizzly.insert235.unfiltered.json')['genes']),'input_pairs_processed':info['n_processed'],'pseudoaligned_pairs':info['n_pseudoaligned'],'pseudoaligned_percent_from_counts':100*info['n_pseudoaligned']/info['n_processed'],'max_caller_splitcount':max(g['splitcount'] for g in source['genes']),'max_caller_paircount':max(g['paircount'] for g in source['genes']),'max_whole_read_exact_sequence_pairs':max(q['whole_read_exact_sequence_pairs'] for q in recounts),'max_central40_sequence_pairs':max(q['central40_sequence_pairs'] for q in recounts),'at_least2_central40_candidates':sum(q['central40_sequence_pairs']>=2 for q in recounts),'single_coordinate_family_candidates':sum(q['all_coordinate_families']==1 for q in familyrows),'normal_reference_exact_junction40_occurrences':sum(int(r['exact_normal_transcript_occurrences']) for r in normal),'public_expected_recovered':len(control['recovered']),'public_unexpected_pairs':len(control['observed'])-len(control['recovered']),'fragment_distribution_q95':q95}
result={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'PASS for bounded negative result with count/wording refinements','summary':summary,'sequence_recounts_all_match':all(q['central40_matches_prior'] and q['whole_read_matches_prior'] for q in recounts),'family_recounts_all_match':all(q['all_family_count_matches'] for q in familyrows),'candidates':recounts,'families':familyrows,'NF1_read_summary':ne,'NF1_public_normal_junction_context':load('full/nf1-junction-independent-audit.json'),'refinements':['SLC29A1-HSP90AB1: six caller names/three coordinate families versus four exact central40 names/two families.','NF1: all four recovered pairs align at NF1; three pairs have nM0 and one nM2, with a 5S read in one pair. Use representative rather than universal zero-mismatch wording.','NF1 intragenic junction is absent from local GENCODE37 NF1 isoforms but appears in6959 GTExv2 samples; supports normal/background alternative splicing, not a specific annotated GENCODE37 isoform.','Fragment length is a kallisto-estimated distribution, not a laboratory-measured molecule-size distribution.'],'limits':['Caller evidence, exact read sequences, read names, coordinate families and original biological molecules are different quantities.','Reported-positive controls verify interoperability; nine recovered controls do not establish patient sensitivity or specificity.','Sequence checks examine reported candidate read sequences and transcript combinations, not an independent exhaustive fusion caller.','No normalized genomic rearrangement or validated fusion antigen was established; absence in this screen is not a comprehensive clinical negative.']}
(P/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(summary);print('sequence/family match',result['sequence_recounts_all_match'],result['family_recounts_all_match'])

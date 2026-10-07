import pathlib,json,collections,pysam,hashlib
O=pathlib.Path(__file__).parent;B=pathlib.Path.home()/'.local/share/codex/caris-analysis';P=json.loads((O/'pairs.json').read_text());S=json.loads((O/'scored-results.json').read_text());rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
# Best orientation/model for every distinct mate sequence; preserve mapping along model rather than infer DNA insert length.
byseq=collections.defaultdict(list)
for r in S['results']:
 if r['query']['kind']=='patient':byseq[r['query']['sequence']].append(r)
pairs=[]
for p in P:
 d={k:p[k]for k in('pair_id','sequence_family','caller_named','Q20_50nt_support','hits')};s1=p['R1_sequence'];s2=rc(p['R2_sequence']);d['read_lengths']=[len(s1),len(s2)];d['mates_exact_reverse_complements']=s1==s2;d['mate_sequence_mismatches_same_length']=[]
 if len(s1)==len(s2):
  for i,(a,b)in enumerate(zip(s1,s2)):
   if a!=b:d['mate_sequence_mismatches_same_length'].append({'R1_position0':i,'R1_base':a,'R1_BQ':ord(p['R1_quality'][i])-33,'R2_position0':len(s2)-1-i,'R2_base_reverse_complemented':b,'R2_BQ':ord(p['R2_quality'][len(s2)-1-i])-33})
 maps=[];mapped_sets=[]
 for mate in(1,2):
  options=byseq[p[f'R{mate}_sequence']];x=max(options,key=lambda x:max(m['score']for m in x['proposed_models']));m=max(x['proposed_models'],key=lambda m:m['score']);al=m['alignment'];co=al['coordinates'];targetpos=[]
  for (ts,qs),(te,qe)in zip(zip(co[0],co[1]),zip(co[0][1:],co[1][1:])):
   if te>ts and qe>qs:targetpos.extend(range(ts,te))
  mapped_sets.append(set(targetpos));maps.append({'mate':mate,'query_id':x['query']['query_id'],'orientation_to_proposed_RNA':x['query']['orientation'],'model':m['id'],'score':m['score'],'query_length':len(x['query']['sequence']),'model_aligned_start0':min(targetpos),'model_aligned_end0':max(targetpos)+1,'query_aligned_fraction':al['aligned_query_fraction'],'crosses_join12_each':al['crosses_join_with12bp_each_side'],'mismatch_count':len(al['mismatches']),'terminal_unaligned_bases':al['query_terminal_unaligned_left']+al['query_terminal_unaligned_right']})
 d['mate_proposed_maps']=maps;d['opposite_mate_orientations']=maps[0]['orientation_to_proposed_RNA']!=maps[1]['orientation_to_proposed_RNA'];d['proposed_RNA_fragment_span']=[min(x['model_aligned_start0']for x in maps),max(x['model_aligned_end0']for x in maps)];d['model_interval_overlap_bp']=max(0,min(x['model_aligned_end0']for x in maps)-max(x['model_aligned_start0']for x in maps));d['aligned_reference_base_overlap_bp']=len(mapped_sets[0]&mapped_sets[1]);pairs.append(d)
# All source STAR records, retaining microhomology counts and alternative-coordinate shifts.
nameid={p['name']:p['pair_id']for p in P};star=[]
for line in(O/'selected-STAR-Chimeric.out.junction').open():
 z=line.strip().split('\t')
 if len(z)<14 or not z[1].isdigit():continue
 star.append({'pair_id':nameid.get(z[9]),'donor_chr':z[0],'donor_intronic_pos1':int(z[1]),'donor_strand':z[2],'acceptor_chr':z[3],'acceptor_intronic_pos1':int(z[4]),'acceptor_strand':z[5],'junction_type':int(z[6]),'repeat_left':int(z[7]),'repeat_right':int(z[8]),'segment1_pos1':int(z[10]),'segment1_CIGAR':z[11],'segment2_pos1':int(z[12]),'segment2_CIGAR':z[13],'num_chimeric_alignments':int(z[14])if len(z)>14 else None,'merged_read_length':int(z[15])if len(z)>15 else None,'STAR_scores':z[16:]})
# Verify microhomology by full250nt joined reference window equality as coordinates shift.
f=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'));bp1=69651217;bp2=125013572
# Fixed outer coordinates keep the reconstructed transcript comparable at every shifted split.
outer1=bp1+100;outer2=bp2-100
construct=lambda delta:rc(f.fetch('chr11',bp1-delta-1,outer1).upper())+rc(f.fetch('chr3',outer2,bp2-delta).upper())
base=construct(0);equiv=[d for d in range(-50,51)if construct(d)==base]
# Shared sequence reassigned between arms by shifting split, not suffix-prefix identity heuristic.
if max(equiv)>0:shared=rc(f.fetch('chr3',bp2-max(equiv),bp2).upper())
elif min(equiv)<0:shared=rc(f.fetch('chr3',bp2,bp2-min(equiv)).upper())
else:shared=''
report={'pairs':pairs,'STAR_records':star,'microhomology':{'method':'Reconstruct same fixed-outer-boundary joined reference sequence while shifting both minus-strand breakpoints together; every delta tested-50..50','equivalent_shifted_split_deltas':equiv,'equivalent_representation_count':len(equiv),'shared_sequence_reassigned':shared,'shared_length':max(equiv)-min(equiv),'original_breakpoints':['chr11:69651217','chr3:125013572'],'coordinate_meaning':'positive delta moves both minus-strand coordinates downward; structural breakpoint is not uniquely localized inside shared sequence'},'summary':{'pairs':len(P),'Q20_support_pairs':sum(p['Q20_50nt_support']for p in P),'Q20_pairs_exact_complementary_mates':sum(p['Q20_50nt_support']and p['mates_exact_reverse_complements']for p in pairs),'Q20_sequence_families':len({p['sequence_family']for p in P if p['Q20_50nt_support']}),'Q20_distinct_model_start_end_patterns':len({tuple(p['proposed_RNA_fragment_span'])for p in pairs if p['Q20_50nt_support']}),'all_distinct_model_start_end_patterns':len({tuple(p['proposed_RNA_fragment_span'])for p in pairs}),'all_pairs_opposite_orientation':all(p['opposite_mate_orientations']for p in pairs)},'limits':['Proposed-RNA geometry is conditional on that hypothesis and is not a genomic DNA insert-size measurement.','Completely overlapping mates observe the same insert; different read names, sequence families and start/end patterns are not original independent molecules.','Microhomology can support alternative breakpoint descriptions or template switching; it does not prove library artifact or a biological rearrangement.']}
(O/'pair-geometry-and-microhomology.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'summary':report['summary'],'microhomology':report['microhomology']}))

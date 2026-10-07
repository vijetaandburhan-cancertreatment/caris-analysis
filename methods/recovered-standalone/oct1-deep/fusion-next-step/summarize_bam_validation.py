import collections,csv,json,pathlib
r=pathlib.Path(__file__).resolve().parent/'full'
a=json.load(open(r/'candidate-bam-audit.json'));e=json.load(open(r/'pizzly.insert235.exact-junctions.json'));p=json.load(open(r/'pizzly.insert235.json'))
pe={(g['geneA']['name'],g['geneB']['name']):g for g in p['genes']};ee={(g['geneA'],g['geneB']):g for g in e['candidate_details']}
rows=[];details=[]
for g in a['candidates']:
 key=(g['geneA']['name'],g['geneB']['name']);reads={(x['read1']['name'],mate):x[mate]for x in pe[key]['readpairs']for mate in ['read1','read2']};events={}
 for t in ee[key]['transcripts']:
  for z in t['evidence']:
   name=z['pair_name'][0];mate=int(z['mate'][-1]);seq=reads[(name,z['mate'])]['seq'];records=[b for b in g['bam_evidence'].get(name,[])if b['mate']==mate and b['forward_sequence']==seq]
   k=(name,mate,z['orientation'],z['left_anchor']);quality=[]
   for b in records:
    q=b['forward_base_qualities']
    if q is None:continue
    if z['orientation']=='reverse-complement':q=q[::-1]
    join=z['left_anchor'];quality.append({'mapq':b['mapq'],'cigar':b['cigar'],'minQ_junction40':min(q[join-20:join+20]),'meanQ_junction40':sum(q[join-20:join+20])/40,'primary':not b['secondary']and not b['supplementary']})
   events[k]={'name':name,'mate':mate,'orientation':z['orientation'],'left_anchor':z['left_anchor'],'right_anchor':z['right_anchor'],'original_BAM_same_sequence':bool(records),'original_BAM_QC':quality}
 names=set();seqs=set()
 for ev in events.values():
  if any(q['primary']and q['mapq']>=20 and q['minQ_junction40']>=20 for q in ev['original_BAM_QC']):
   names.add(ev['name']);pair=next(z for z in pe[key]['readpairs']if z['read1']['name']==ev['name']);seqs.add((pair['read1']['seq'],pair['read2']['seq']))
 row={'geneA':key[0],'geneB':key[1],'raw_name_pairs':len(g['names']),'exact20_name_pairs':ee[key]['exact20nt_join_name_pairs'],'exact20_sequence_pairs':ee[key]['exact20nt_join_sequence_pairs'],'BAM_verified_exact_Q20_MAPQ20_name_pairs':len(names),'BAM_verified_exact_Q20_MAPQ20_sequence_pairs':len(seqs)};rows.append(row);details.append(dict(row,read_mate_checks=list(events.values()),junctions=g['junction_annotation']))
res={'interpretation':'Direct evidence review of weak candidate reads. Passing Q/MAPQ checks does not validate a biological or somatic fusion.','caveats':['BAM mapping quality applies to the single aligned portion; it does not verify the other fusion partner.','Counts collapse names and identical paired sequences, not optical/PCR molecules or cell clones.','Reference-transcript/junction checks cannot replace a full genome aligner and orthogonal junction validation.'],'summary':rows,'details':details}
(r/'candidate-bam-validation.json').write_text(json.dumps(res,indent=2)+'\n')
with open(r/'candidate-bam-validation.tsv','w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))

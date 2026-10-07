"""Coordinate-family check for candidate read pairs; no UMI molecule claim."""
import collections,json,pathlib,re
r=pathlib.Path(__file__).resolve().parent/'full';d=json.loads((r/'candidate-bam-audit.json').read_text());rows=[];details=[]
w=json.loads((r/'pizzly.insert235.junction-window-validation.json').read_text());window_names={(g['geneA'],g['geneB']):{e['name'] for t in g['junctions'] for e in t['read_evidence']} for g in w['details']}
def prime5(b):
 cig=[(int(n),op)for n,op in re.findall(r'(\d+)([MIDNSHP=X])',b['cigar'])]
 if b['flag']&16:
  clip=0
  for n,op in reversed(cig):
   if op in 'SH':clip+=n
   else:break
  return b['reference'],b['end1']+clip,'reverse'
 clip=0
 for n,op in cig:
  if op in 'SH':clip+=n
  else:break
 return b['reference'],b['start1']-clip,'forward'
for g in d['candidates']:
 families=collections.defaultdict(list);incomplete=[]
 for name,ev in g['bam_evidence'].items():
  prim=[b for b in ev if not b['secondary']and not b['supplementary']and not b['qc_failed']]
  m1=[b for b in prim if b['mate']==1];m2=[b for b in prim if b['mate']==2]
  if len(m1)!=1 or len(m2)!=1:incomplete.append(name);continue
  key=(prime5(m1[0]),prime5(m2[0]));families[key].append(name)
 central_names=window_names.get((g['geneA']['name'],g['geneB']['name']),set())
 central_families=sum(bool(set(v)&central_names) for v in families.values())
 row={'central40_supported_names':len(central_names),'central40_primary_pair_coordinate_families':central_families,'geneA':g['geneA']['name'],'geneB':g['geneB']['name'],'candidate_names':len(g['names']),'primary_pair_unclipped5prime_coordinate_families':len(families),'names_without_complete_primary_pair':len(incomplete)}
 rows.append(row);details.append(dict(row,families=[{'mate1_chrom_unclipped5prime_strand':k[0],'mate2_chrom_unclipped5prime_strand':k[1],'names':v}for k,v in families.items()],incomplete=incomplete))
result={'method':'Group names by both mates chromosome, strand and unclipped5prime genomic coordinate, as read from original STAR BAM. This uses known alignment portions; it is not a UMI count or re-alignment.','interpretation':'Distinct names or sequence pairs can still be PCR/optical families. A single coordinate family is insufficient independent-template support for a new clinical fusion call. Different coordinates do not prove biological molecules or a rearrangement.','summary':rows,'details':details}
(r/'candidate-fragment-families.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(rows,indent=2))

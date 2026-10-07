#!/usr/bin/env python3
"""Recover explicit breakend-proximity STAR names and conservative reference-junction hypotheses.
This does not recover Arriba's filtered read list, nor validate a DNA rearrangement.
"""
import collections,copy,json,pathlib,re,pysam
base=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-fusion-read-audit');out=base/'supplemental-inventory';out.mkdir(exist_ok=True)
rows=json.load(open(base/'patient-inventory/review-candidates.json')); genome='/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference/GRCh38.primary_assembly.genome.fa'; fa=pysam.FastaFile(genome)
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def endpoint(bp,d):
 chrom,pos=bp.split(':');return chrom,d,int(pos)
def ordered(a,b):return sorted((a,b),key=lambda x:(x[0],x[1],x[2]))
def seq(chrom,d,pos,side,arm):
 s=fa.fetch(chrom,pos-arm,pos) if d=='downstream' else fa.fetch(chrom,pos-1,pos+arm-1)
 if (side=='left' and d=='upstream') or (side=='right' and d=='downstream'):s=rc(s)
 return s.upper()
sup=[];verify=[]
for row in rows:
 ep=[endpoint(row[f'breakpoint{i}'],row[f'direction{i}']) for i in [1,2]]
 patterns=[{'arm_bases':a,'sequence':seq(*ep[0],'left',a)+seq(*ep[1],'right',a)} for a in [15,20,25]]
 if row['junction_patterns']:verify.append({'id':row['candidate_id'],'gene1':row['gene1'],'gene2':row['gene2'],'reference_exactly_matches_caller_at_arms':[a['arm_bases'] for a in patterns if any(a==b for b in row['junction_patterns'])]})
 if row['source_set']=='discarded' or row['gene1']=='MADCAM1':
  x=copy.deepcopy(row);x['original_source_set']=row['source_set'];x['source_set']='reference_reconstructed_'+row['source_set'];x['original_distinct_named_reads']=row['distinct_named_reads'];x['original_read_names']=row['read_names'];x['junction_patterns']=patterns;x['junction_hypothesis']='Exact GRCh38 genomic sequence joined at reported breakpoints, no splice/SNV/insertion reconstruction. Negative exact counts are uninformative for alternate RNA sequence.';x['read_names']=[];x['name_source']='STAR chimeric lines with same two genomic breakend directions, each within25bp of reported breakpoints; not Arriba filtered names';x['_endpoints']=ordered(*ep);x['_star_matches']=[];sup.append(x)
 elif row['fusion_transcript'].count('|')==2:
  left,insert,right=row['fusion_transcript'].split('|');l=re.search('[ACGTacgt]+$',left);r=re.match('[ACGTacgt]+',right)
  if l and r and re.fullmatch('[ACGTacgt]+',insert):
   x=copy.deepcopy(row);x['source_set']='insertion_consensus_'+row['source_set'];x['original_source_set']=row['source_set'];x['junction_patterns']=[{'arm_bases':a,'sequence':(l[0][-a:]+insert+r[0][:a]).upper(),'inserted_bases':insert.upper()} for a in [15,20,25] if min(len(l[0]),len(r[0]))>=a];x['junction_hypothesis']='Exact caller-assembled left arm, full non-template insertion, right arm; all marker bases quality-filtered.';sup.append(x)
groups=collections.defaultdict(list)
for r in sup:
 if '_endpoints' in r:
  e=r['_endpoints'];groups[tuple((x[0],x[1]) for x in e)].append(r)
star='work/oct4-followup/genome-fusion/patient-D8/STAR.Chimeric.out.junction';n=0
with open(star) as f:
 for line in f:
  z=line.rstrip().split('\t')
  if len(z)<14 or not z[1].isdigit():continue
  donor=(z[0],'downstream' if z[2]=='+' else 'upstream',int(z[1])-1 if z[2]=='+' else int(z[1])+1)
  acceptor=(z[3],'upstream' if z[5]=='+' else 'downstream',int(z[4])+1 if z[5]=='+' else int(z[4])-1)
  ep=ordered(donor,acceptor);key=tuple((x[0],x[1]) for x in ep)
  for row in groups.get(key,[]):
   shifts=[observed[2]-expected[2] for observed,expected in zip(ep,row['_endpoints'])]
   if all(abs(d)<=25 for d in shifts):
    row['_star_matches'].append({'name':z[9],'shifts':shifts,'exact_breakend_coordinates':all(d==0 for d in shifts),'star_columns':z});n+=1
for row in sup:
 if '_star_matches' in row:
  row['read_names']=sorted(set(m['name'] for m in row['_star_matches']));row['distinct_named_reads']=len(row['read_names']);row['star_exact_coordinate_names']=sorted(set(m['name'] for m in row['_star_matches'] if m['exact_breakend_coordinates']))
(out/'review-candidates.json').write_text(json.dumps(sup,indent=2)+'\n')
(out/'coordinate-formula-accepted-controls.json').write_text(json.dumps(verify,indent=2)+'\n')
print(json.dumps({'candidates':len(sup),'by_type':dict(collections.Counter(x['source_set'] for x in sup)),'star_matching_lines':n,'derived_unique_names':len(set(n for r in sup for n in r['read_names'])),'accepted_reference_controls':verify},indent=2))

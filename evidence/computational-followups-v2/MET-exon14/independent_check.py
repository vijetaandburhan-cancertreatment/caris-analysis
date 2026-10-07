#!/usr/bin/env python3
"""Independent MET exon14 annotation/CIGAR check; no production code imports."""
import pathlib,json,re,hashlib,pysam,collections,copy,shutil
O=pathlib.Path(__file__).resolve().parent;B=pathlib.Path.home()/'.local/share/codex/caris-analysis';D=B/'oct5-MET-exon14-audit';D.mkdir(exist_ok=True)
G=B/'genome-fusion-reference/gencode.v37.primary_assembly.annotation.gtf';fa=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'));tx='ENST00000397752.8';ex={};cds=[];txline=''
for line in G.open():
 if line.startswith('#')or tx not in line:continue
 c=line.rstrip().split('\t');attrs=dict(re.findall(r'(\w+) "([^"]*)"',c[8]))
 if attrs.get('transcript_id')!=tx:continue
 if c[2]=='transcript':txline=line.strip()
 if c[2]=='exon':ex[int(re.search(r'exon_number \"?(\d+)',c[8]).group(1))]={'chrom':c[0],'start0':int(c[3])-1,'end0':int(c[4]),'strand':c[6]}
 if c[2]=='CDS':cds.append((int(c[3])-1,int(c[4])))
assert len(ex)==21 and all(x['strand']=='+'for x in ex.values())and 'MANE_Select'in txline
assert all(ex[i]['end0']<ex[i+1]['start0']for i in range(1,21))
junctions={'exon13_to14':(ex[13]['end0'],ex[14]['start0']),'exon14_to15':(ex[14]['end0'],ex[15]['start0']),'exon13_to15_skip14':(ex[13]['end0'],ex[15]['start0'])}

def permitted(r):
 return not(r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail)and r.mapping_quality>=20 and(not r.has_tag('NH')or r.get_tag('NH')==1)
def walker(r):
 if not permitted(r)or not r.query_sequence or r.query_qualities is None:return set()
 rp=r.reference_start;qp=0;mapping={};skips=[];cg=r.cigartuples or[]
 for i,(op,n)in enumerate(cg):
  if op in(0,7,8):mapping.update((rp+k,qp+k)for k in range(n))
  if op==3:skips.append((rp,rp+n,i,qp))
  if op in(0,2,3,7,8):rp+=n
  if op in(0,1,4,7,8):qp+=n
 out=set()
 for s,e,i,q in skips:
  if i==0 or i+1==len(cg)or cg[i-1][0]not in(0,7,8)or cg[i+1][0]not in(0,7,8)or min(cg[i-1][1],cg[i+1][1])<12:continue
  pos=list(range(s-12,s))+list(range(e,e+12));qq=[mapping.get(p)for p in pos]
  if qq!=list(range(q-12,q+12)):continue
  if min(r.query_qualities[k]for k in qq)<20:continue
  obs=''.join(r.query_sequence[k]for k in qq);ref=fa.fetch('chr7',s-12,s).upper()+fa.fetch('chr7',e,e+12).upper()
  if obs==ref:out.add((s,e))
 return out
header=pysam.AlignmentHeader.from_dict({'SQ':[{'SN':'chr7','LN':159345973}]});s,e=junctions['exon13_to15_skip14']
def synth(label,alen=50,blen=50,op=3,flag=0,mapq=60,nh=1,qual=35,mismatch=False):
 r=pysam.AlignedSegment(header);r.query_name=label;r.reference_id=0;r.reference_start=s-alen;r.flag=flag;r.mapping_quality=mapq;r.cigartuples=[(0,alen),(op,e-s),(0,blen)];seq=fa.fetch('chr7',s-alen,s).upper()+fa.fetch('chr7',e,e+blen).upper()
 if mismatch:seq=seq[:alen-1]+('C'if seq[alen-1]!='C'else'A')+seq[alen:]
 r.query_sequence=seq;r.query_qualities=[qual]*len(seq);r.set_tag('NH',nh);return r
cases=[('positive',{},True),('reverse_stored_orientation',{'flag':16},True),('11bp_anchor',{'alen':11},False),('anchor_mismatch',{'mismatch':True},False),('low_BQ',{'qual':19},False),('deletion_not_splice',{'op':2},False),('secondary',{'flag':256},False),('multi_mapper',{'nh':2},False),('MAPQ19',{'mapq':19},False)]
controls=[]
for name,kwargs,want in cases:
 r=synth(name,**kwargs);got=(s,e)in walker(r);assert got==want;controls.append({'case':name,'expected_skip':want,'observed_skip':got,'CIGAR':r.cigarstring})
counts={k:set()for k in junctions};rgcounts={k:set()for k in junctions};all_js=collections.Counter();nall=npass=0
bam=pysam.AlignmentFile(str(B/'TN26-279853/RNA_TN26-279853.bam'),'rb',index_filename=str(pathlib.Path.cwd()/'work/oct1-analysis/RNA_TN26-279853.bam.bai'))
for r in bam.fetch('chr7',ex[13]['start0']-200,ex[15]['end0']+200):
 nall+=1;npass+=permitted(r);js=walker(r)
 for k,j in junctions.items():
  if j in js:counts[k].add(r.query_name);rgcounts[k].add((r.query_name,r.get_tag('RG')if r.has_tag('RG')else''))
 for j in js:all_js[j]+=1
prod=json.loads((O/'junction-check.json').read_text());sjrows={};sjfile=O.parent/'patient-D8/STAR.SJ.out.tab'
for line in sjfile.open():
 a=line.split()
 if a[0]!='chr7':continue
 start0=int(a[1])-1;end0=int(a[2]);j=(start0,end0)
 for k,target in junctions.items():
  if j==target:sjrows[k]={'STAR_start1':int(a[1]),'STAR_end1':int(a[2]),'converted_interval0':[start0,end0],'strand':int(a[3]),'motif':int(a[4]),'annotated':int(a[5]),'unique_reads':int(a[6]),'multi_reads':int(a[7]),'max_overhang':int(a[8])}
result={'status':'pass','transcript':tx,'MANE_Select_verified_in_G37':True,'exon_count':len(ex),'exons13_14_15':{k:ex[k]for k in(13,14,15)},'exon14_length_bp':ex[14]['end0']-ex[14]['start0'],'exon14_CDS_intersection_bp':sum(max(0,min(b,ex[14]['end0'])-max(a,ex[14]['start0']))for a,b in cds),'old_BAM_alignments_fetched':nall,'old_BAM_filter_passing_alignments':npass,'results':{},'synthetic_controls':controls,'public_source':'https://pubmed.ncbi.nlm.nih.gov/28522754/','limits':['Synthetic controls validate this counter and coordinates, not clinical assay sensitivity.','New STAR.SJ counts are aligner read counts, not directly comparable to the old strict query-name counts.','No detected canonical skip junction does not exclude every MET alteration or a low-abundance/variant-junction transcript.','No treatment recommendation or MET inhibitor eligibility is inferred.']}
for k,j in junctions.items():
 p=prod['results'][k];assert set(p['original_strict_names'])==counts[k];assert(j[0],j[1])==(p['intron_start0'],p['intron_end0'])
 marker=fa.fetch('chr7',j[0]-25,j[0]).upper()+fa.fetch('chr7',j[1],j[1]+25).upper();assert marker==p['RNA_junction_25bp_each_side']
 result['results'][k]={'interval0':j,'old_strict_query_names':len(counts[k]),'old_strict_query_name_RG_fragments':len(rgcounts[k]),'exact_query_name_set_matches_production':True,'new_SJ':sjrows.get(k),'independent_marker50':marker,'splice_motif':fa.fetch('chr7',j[0],j[0]+2).upper()+'-'+fa.fetch('chr7',j[1]-2,j[1]).upper()}
for x in(O,D):(x/'independent-junction-check.json').write_text(json.dumps(result,indent=2)+'\n')
shutil.copy2(__file__,D/pathlib.Path(__file__).name)
print(json.dumps(result))

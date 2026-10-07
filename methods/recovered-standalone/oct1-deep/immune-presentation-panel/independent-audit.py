"""Read-only arithmetic/provenance audit with bounded coding-depth checks."""
from pathlib import Path
from collections import defaultdict,Counter
import json,csv,gzip,re,datetime
import pysam
P=Path(__file__).resolve().parent;ROOT=P.parents[2];RAW=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
A=json.loads((P/'target-annotation.json').read_text())['genes'];SCREEN=json.loads((P/'screen-results.json').read_text());ALLANN={r['transcript']:r for r in json.loads((ROOT/'work/oct1-deep/coding-indels/selected-coding-transcripts.json').read_text())};PUB=json.loads((P/'extended-protein-candidate-public-annotation.json').read_text())['results']
BASES='TCAG';AA='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';CODE=dict(zip((a+b+c for a in BASES for b in BASES for c in BASES),AA));COMP=str.maketrans('ACGT','TGCA')
def translation(s):return ''.join(CODE[s[i:i+3]] for i in range(0,len(s)-2,3))
def merge(intervals):
 pts=[]
 for s,e in sorted(intervals):
  if pts and s<=pts[-1][1]:pts[-1][1]=max(e,pts[-1][1])
  else:pts.append([s,e])
 return pts
annrows=[];cons=[];pubrows=[]
for gene,a in A.items():
 ann=ALLANN[a['selected_transcript']];annrows.append({'gene':gene,'transcript':a['selected_transcript'],'is_MANE':'MANE_Select' in ann['tags'],'priority':ann['priority'],'CDS_matches_source':sorted(a['selected_CDS'])==sorted(ann['CDS'])})
 w=json.loads((P/f'{gene}.hg38-reference.json').read_text());ref=w['dna'].upper();gpos=[p for s,e in sorted(a['selected_CDS']) for p in range(s,e)];seq=''.join(ref[p-w['start']] for p in gpos)
 if a['strand']=='-':seq=seq.translate(COMP)[::-1];gpos.reverse()
 protein=translation(seq)
 for v in SCREEN['candidates']:
  if v['gene']!=gene:continue
  assert len(v['ref'])==len(v['alt'])==1
  idx=gpos.index(v['pos1']-1);refbase=v['ref'] if a['strand']=='+' else v['ref'].translate(COMP);altbase=v['alt'] if a['strand']=='+' else v['alt'].translate(COMP);assert seq[idx]==refbase
  mut=seq[:idx]+altbase+seq[idx+1:];mp=translation(mut);aa=idx//3;pred='=' if mp[aa]==protein[aa] else protein[aa]+str(aa+1)+mp[aa]
  cons.append({'gene':gene,'pos1':v['pos1'],'protein':pred,'matches':pred==v['protein'],'reference_no_internal_stop':protein.startswith('M') and '*' not in protein})
for rr in PUB:
 r=rr['candidate'];v=next(v for v in SCREEN['candidates'] if v['gene']==r['gene'] and v['pos1']==r['position_1based']);source='Caris VCF' if v['Caris_VCF_exact_matches'] else 'research read-count screen';expected=float(v['Caris_VCF_exact_matches'][0]['VAF']) if v['Caris_VCF_exact_matches'] else v['screen_fraction'];matches=[];freq=[]
 for ident in rr['variation_records_cached']:
  d=json.loads((P/'public-population'/f'{ident}.variation.json').read_text())['result'];m=[m for m in d['mappings'] if m.get('assembly_name')=='GRCh38' and m.get('seq_region_name')==r['chrom'][3:] and m['start']==m['end']==r['position_1based'] and m['strand']==1 and r['ref'] in m['allele_string'].split('/') and r['alt'] in m['allele_string'].split('/')];matches.extend(m)
  freq.extend(q for q in d.get('populations',[]) if q.get('allele')==r['alt'] and q.get('population') in ['gnomADe:ALL','gnomADg:ALL','1000GENOMES:phase_3:ALL'] and q.get('frequency') is not None)
 pubrows.append({'gene':r['gene'],'pos1':r['position_1based'],'public_mapping_pass':bool(matches),'global_frequency_at_least_1percent':max(float(q['frequency']) for q in freq)>=.01,'DNA_fraction_source_pass':r['DNA_allele_fraction_source']==source and abs(float(r['DNA_allele_fraction'])-expected)<1e-12})
# Recount one small set of gene intervals and candidate positions independently;
# retains exactly the source read-depth definition, not molecule depth.
coverage=[]
with pysam.AlignmentFile(str(RAW/'DNA_TN26-279853.bam'),'rb',index_filename=str(ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai')) as bam:
 for gene in ['B2M','STAT1','RFXAP','NLRC5']:
  a=A[gene];regions=merge([(s-4,e+4) for s,e in a['selected_CDS']]);depth={}
  for s,e in regions:
   arr=bam.count_coverage(a['chrom'],s,e,quality_threshold=20,read_callback=lambda r:not r.flag&(4|256|512|1024|2048) and r.mapping_quality>=20)
   for j,q in enumerate(zip(*arr)):depth[s+j]=sum(q)
  old=next(r for r in SCREEN['gene_summaries'] if r['gene']==gene);new={'CDS_plus4_bases':len(depth),'bases_ge20':sum(d>=20 for d in depth.values()),'minimum_strict_read_depth':min(depth.values())};low=[{'pos1':p+1,'read_depth':d,'inside_selected_CDS':any(s<=p<e for s,e in a['selected_CDS'])} for p,d in depth.items() if d<20]
  coverage.append({'gene':gene,**new,'matches':all(new[k]==old[k] for k in new),'low_depth_positions':low})
# Independently sum the expression export by gene, including duplicate symbols.
expr=defaultdict(lambda:[0.,0.])
with (RAW/'RNA_TN26-279853.geneTPM_nodup.csv').open() as f:
 for r in csv.DictReader(f):
  if r['Gene'] in A:expr[r['Gene']][0]+=float(r['TPM']);expr[r['Gene']][1]+=float(r['NumReads'])
with (P/'selected-expression.tsv').open() as f:expression=[{'gene':r['Gene'],'matches':abs(float(r['TPM'])-expr[r['Gene']][0])<1e-8 and abs(float(r['NumReads'])-expr[r['Gene']][1])<1e-8} for r in csv.DictReader(f,delimiter='\t')]
result={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'PASS after source-label corrections','transcripts':annrows,'consequences':cons,'public_allele_audit':pubrows,'representative_depth_audit':coverage,'expression_audit':expression,'summary_recomputed':{'genes':len(A),'MANE':sum(r['is_MANE'] for r in annrows),'candidates':len(SCREEN['candidates']),'types':dict(Counter(r['consequence'] for r in SCREEN['candidates'])),'total_selected_CDS_plus4':sum(r['CDS_plus4_bases'] for r in SCREEN['gene_summaries']),'total_bases_ge20':sum(r['bases_ge20'] for r in SCREEN['gene_summaries'])},'limits':['Read counts do not collapse overlapping mates; they must remain labelled read depth and screen allele fraction.','Initial filters do not establish a validated sensitivity, clonal limit, somatic classification or intact protein function.','HLA-B expression/retention, complete isoform coverage, structural/regulatory/epigenetic changes and later acquired resistance are outside this screen.','Common population alleles were not assigned clinical benignity or proof of germline status.']}
(P/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(result['summary_recomputed']);print('coverage',coverage);print('all consequences/pop/expression pass',all(r['matches'] for r in cons),all(r['public_mapping_pass'] and r['global_frequency_at_least_1percent'] and r['DNA_fraction_source_pass'] for r in pubrows),all(r['matches'] for r in expression))

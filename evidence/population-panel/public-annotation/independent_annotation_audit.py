from pathlib import Path
import csv,json,hashlib,datetime,collections,re
ROOT=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct5-population-panel-concordance-v1/public-annotation')
OUT=Path(__file__).resolve().parent
PANEL=ROOT/'public-panel-gene-annotations.tsv'
GTF=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference/gencode.v37.primary_assembly.annotation.gtf')
sha=lambda p:hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
rows=list(csv.DictReader(PANEL.open(),delimiter='\t'))
groups=collections.defaultdict(list)
for r in rows: groups[(r['chrom'],r['gene_overlap_category'])].append(r)
selected={}
for key,rr in sorted(groups.items()):
 for idx in (0,len(rr)//2,len(rr)-1):
  r=rr[idx];selected[(r['chrom'],int(r['pos1']))]=r
# Add boundaries selected solely from the public annotation inventory, with adjacent
# panel positions immediately below/above each named start/end, independent of alleles.
inv=list(csv.DictReader((ROOT/'overlapping-gene-inventory.tsv').open(),delimiter='\t'))
names={'BAP1','RASA1','MTAP','CDKN2A','CDKN2B','ACTB','TP53','LATS2','MYC','RB1'}
for g in inv:
 if g['gene_symbol'] not in names:continue
 rr=[r for r in rows if r['chrom']==g['chrom']]
 for b in (int(g['start1']),int(g['end1'])):
  for direction in (-1,1):
   choices=[r for r in rr if (int(r['pos1'])-b)*direction>=0]
   if choices:
    r=min(choices,key=lambda r:abs(int(r['pos1'])-b));selected[(r['chrom'],int(r['pos1']))]=r
expected={k:{'genes':set(),'pc_exons':set(),'pct_exons':set()} for k in selected}
byc=collections.defaultdict(list)
for c,p in selected:byc[c].append(p)
meta={};h=hashlib.sha256()
with GTF.open('rb') as f:
 for line in f:
  h.update(line)
  if line.startswith(b'#'):continue
  z=line.decode().rstrip('\n').split('\t')
  if z[0] not in byc or z[2] not in ('gene','exon'):continue
  start,end=int(z[3]),int(z[4]);at=dict(re.findall(r'(\w+) "([^"]*)"',z[8]));gid=at['gene_id']
  if z[2]=='gene':meta[gid]=at
  for p in byc[z[0]]:
   if not (start<=p<=end):continue
   e=expected[(z[0],p)]
   if z[2]=='gene':e['genes'].add(gid)
   elif at['gene_type']=='protein_coding':
    e['pc_exons'].add(gid)
    if at.get('transcript_type')=='protein_coding':e['pct_exons'].add(gid)
def decode(x):return set() if x=='.' else set(x.split(';'))
dif=[];evidence=[]
for k,r in sorted(selected.items()):
 e=expected[k]; checks={
  'overlap_gene_ids':e['genes']==decode(r['overlap_gene_ids']),
  'pc_gene_exon_gene_ids':e['pc_exons']==decode(r['pc_gene_exon_gene_ids']),
  'pc_transcript_exon_gene_ids':e['pct_exons']==decode(r['pc_transcript_exon_gene_ids']),
  'unique_gene_id':r['unique_gene_id']==(next(iter(e['genes'])) if len(e['genes'])==1 else '.'),
  'gene_overlap_category':r['gene_overlap_category']==('intergenic' if not e['genes'] else 'unique_gene' if len(e['genes'])==1 else 'multi_gene')}
 evidence.append({'chrom':k[0],'pos1':k[1],'expected_gene_ids':sorted(e['genes']),'expected_pc_exon_gene_ids':sorted(e['pc_exons']),'expected_pc_transcript_exon_gene_ids':sorted(e['pct_exons']),'checks':checks})
 if not all(checks.values()):dif.append(evidence[-1])
result={'status':'PASS' if not dif else 'FAIL','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Public annotation only; no patient BAM or genotype data read.','method':'Independent sequential GTF attribute parsing and direct inclusive containment tests for deterministic chromosome/category cases and public named-gene start/end neighboring panel sites. No primary bisect routine imported.','sampled_loci':len(selected),'annotation_sha256':sha(PANEL),'gtf_sha256':h.hexdigest(),'sampled_category_counts':dict(collections.Counter(r['gene_overlap_category'] for r in selected.values())),'all_locus_category_recount':dict(collections.Counter(r['gene_overlap_category'] for r in rows)),'differences':dif,'checks':evidence,'limitations':'A bounded annotation spot-check, not independent reannotation of all 165782 loci. Gene spans/exons include UTR; labels do not establish patient expression or independence.'}
(OUT/'independent-annotation-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('checks',)},indent=2))

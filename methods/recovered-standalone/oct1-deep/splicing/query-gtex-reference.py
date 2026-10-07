"""Public gene-only queries: no patient IDs, sequences, genotypes, counts or files transmitted."""
from pathlib import Path
import csv,json,urllib.request,urllib.parse,time,concurrent.futures,hashlib,io,collections
OUT=Path(__file__).resolve().parent;D=OUT/'gtex-public';D.mkdir(exist_ok=True)
rows=list(csv.DictReader((OUT/'candidate-reference-audit.tsv').open(),delimiter='\t'))
genes=sorted({r['gene'] for r in rows})
fields='snaptron_id,chromosome,start,end,strand,samples_count,coverage_sum,coverage_avg,coverage_median,annotated'

def fetch(gene):
 params={'regions':gene,'fields':fields,'header':'1'}
 url='https://snaptron.cs.jhu.edu/gtexv2/snaptron?'+urllib.parse.urlencode(params)
 path=D/(gene+'.gtexv2.all.tsv');meta_path=D/(gene+'.gtexv2.all.query.json')
 if path.exists() and meta_path.exists():return gene,list(csv.DictReader(path.open(),delimiter='\t')),json.load(meta_path.open())
 for attempt in range(2):
  try:
   t=time.monotonic();req=urllib.request.Request(url,headers={'User-Agent':'Python public reference gene query'})
   with urllib.request.urlopen(req,timeout=60) as response:
    data=response.read(8_000_001);assert len(data)<=8_000_000,'Unexpected response >8MB; stop instead of fetching larger data'
   text=data.decode();table=list(csv.DictReader(io.StringIO(text),delimiter='\t'));assert 'samples_count' in table[0] if table else 'DataSource:Type' in text
   path.write_bytes(data);meta={'gene':gene,'compilation':'gtexv2','url':url,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'returned_junction_rows':len(table),'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'elapsed_seconds':time.monotonic()-t,'transmitted':'public gene symbol and output field names only; no patient sequences/files/identifiers/variant list'}
   meta_path.write_text(json.dumps(meta,indent=2)+'\n');print(gene,len(table),'rows',len(data),'bytes',flush=True)
   return gene,table,meta
  except Exception as exc:
   if attempt==1:return gene,[],{'gene':gene,'url':url,'error':repr(exc)}
   time.sleep(2)
raw={}; metas={}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for gene,table,meta in pool.map(fetch,genes):raw[gene]=table;metas[gene]=meta
index={gene:{(r['chromosome'],int(r['start'])-1,int(r['end']),r['strand']):r for r in table} for gene,table in raw.items()}
result=[]
for r in rows:
 k=r['chrom'],int(r['intron_start0']),int(r['intron_end0']),r['target_strand']; hit=index[r['gene']].get(k)
 out={**r,'GTEx_query_success':int('error' not in metas[r['gene']]),'GTEx_coordinate_strand_match':int(hit is not None),'GTEx_samples_with_junction':hit['samples_count'] if hit else '', 'GTEx_junction_read_sum':hit['coverage_sum'] if hit else '', 'GTEx_mean_reads_in_positive_samples':hit['coverage_avg'] if hit else '', 'GTEx_median_reads_in_positive_samples':hit['coverage_median'] if hit else '', 'GTEx_snaptron_id':hit['snaptron_id'] if hit else '', 'GTEx_annotation_labels':hit['annotated'] if hit else ''}
 result.append(out)
with (OUT/'candidate-gtex-comparison.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
summary={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'query_genes':len(genes),'successful_gene_queries':sum('error' not in m for m in metas.values()),'total_response_bytes':sum(m.get('bytes',0) for m in metas.values()),'candidates':len(result),'exact_coordinate_strand_matches':sum(r['GTEx_coordinate_strand_match'] for r in result),'reference_aware_passing_candidates':sum(r['reference_aware_research_pass']=='1' for r in result),'reference_aware_passing_with_GTEx_match':sum(r['reference_aware_research_pass']=='1' and r['GTEx_coordinate_strand_match'] for r in result),'queries':metas,'coordinate_conversion':'Snaptron hg38 intron coordinates are 1-based inclusive: compare start-1 and end to local 0-based half-open CIGAR N intervals, with same strand.','limitations':['GTEx contains public normal-tissue sequencing, not a matched normal from this patient.','Observed junction alignments can reflect normal isoforms, background mis-splicing, or shared technical artifacts. Their presence argues against tumor-specificity but does not prove a particular biological mechanism.','Absence from this returned atlas is not proof of tumor-specificity; tissue representation, expression, coverage and processing differ.','Counts are GTEx read-based aggregate statistics; not directly comparable to this sample fragment counts and not calibrated transcript PSI.','No full per-sample genotype/metadata or protein/T-cell measurements were retrieved.'],'documentation':['https://snaptron.idies.jhu.edu/data.html','https://snaptron.idies.jhu.edu/reftables.html','https://snaptron.idies.jhu.edu/wsi.html']}
(OUT/'gtex-comparison-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='queries'},indent=2),flush=True)

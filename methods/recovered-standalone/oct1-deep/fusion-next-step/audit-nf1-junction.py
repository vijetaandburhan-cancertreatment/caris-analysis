"""Check NF1 intragenic junction against local annotation and public GTEx only."""
from pathlib import Path
from collections import defaultdict
import gzip,re,json,urllib.request,datetime,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parents[2];start,end=31343135,31345388;exons=defaultdict(list);meta={}
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  a=line.rstrip().split('\t')
  if a[0]!='chr17' or a[2]!='exon' or 'gene_name "NF1"' not in a[8]:continue
  tid=re.search('transcript_id "([^"]+)"',a[8])[1];exons[tid].append((int(a[3])-1,int(a[4])));meta[tid]=dict(re.findall('(\w+) "([^"]+)"',a[8]))
hits=[]
for tid,ex in exons.items():
 ex.sort()
 if any(a[1]==start and b[0]==end for a,b in zip(ex,ex[1:])):hits.append({'transcript':tid,'type':meta[tid].get('transcript_type'),'name':meta[tid].get('transcript_name')})
url='https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=NF1&fields=snaptron_id,chromosome,start,end,strand,samples_count,coverage_sum,coverage_avg,coverage_median,annotated&header=1';cache=P/'full/NF1.gtexv2.audit.tsv';rows=[];error=None
try:
 if not cache.exists():
  with urllib.request.urlopen(url,timeout=45) as r:data=r.read()
  cache.write_bytes(data)
 lines=cache.read_text().splitlines();header=next(x for x in lines if 'snaptron_id' in x and 'chromosome' in x).split('\t')
 for line in lines:
  a=line.split('\t')
  if len(a)!=len(header) or a==header:continue
  row=dict(zip(header,a))
  try:
   if row['chromosome']=='chr17' and int(row['start'])-1==start and int(row['end'])==end and row['strand']=='+':rows.append(row)
  except (KeyError,ValueError):pass
except Exception as e:error=repr(e)
out={'junction_GRCh38_0based_halfopen':['chr17',start,end,'+'],'GENCODE37_exact_transcript_junction_matches':hits,'gtex_source_url':url,'gtex_rows':rows,'gtex_error':error,'gtex_cache_sha256':hashlib.sha256(cache.read_bytes()).hexdigest() if cache.exists() else None,'interpretation':'An intragenic genomic alignment supplies a competing explanation for the proposed NF1-AK4 fusion. Public normal-junction presence, if found, supports a normal/background splice explanation but does not establish every current read source or tumor specificity. Absence from annotation or public data is not proof of a novel tumor event.'}
(P/'full/nf1-junction-independent-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

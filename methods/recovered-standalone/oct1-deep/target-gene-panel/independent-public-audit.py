import json,pathlib,urllib.request,concurrent.futures,datetime,collections,time
B=pathlib.Path(__file__).resolve().parent; cache=B/'independent-public-cache';cache.mkdir(exist_ok=True)
C=json.loads((B/'screen-results.json').read_text())['candidates']; V=json.loads((B.parents[1]/'oct1-analysis/variants.vcf-records.json').read_text())
index=collections.defaultdict(list)
for v in V:index[(v['CHROM'],int(v['POS']),v['REF'],v['ALT'])].append(v)
all_exact=[]
for c in C:
 m=index[(c['chrom'],c['pos1'],c['ref'],c['alt'])]
 all_exact.append({'gene':c['gene'],'chrom':c['chrom'],'pos1':c['pos1'],'ref':c['ref'],'alt':c['alt'],'exact_matches':len(m),'filters':[v['FILTER'] for v in m],'clinical':[v.get('info',{}).get('CI') for v in m]})
assert all(x['exact_matches'] for x in all_exact)
def get(url,name):
 p=cache/name
 if p.exists():return json.loads(p.read_text())['result']
 req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'Reference-variant-audit/1.0'})
 with urllib.request.urlopen(req,timeout=20) as r:j=json.load(r)
 p.write_text(json.dumps({'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'result':j},indent=2)+'\n');return j
def one(c):
 ch=c['chrom'].replace('chr','');pos=c['pos1'];ref=c['ref'];alt=c['alt'];r={k:c[k] for k in ['gene','chrom','pos1','ref','alt','protein','transcript']}
 try:
  url=f'https://rest.ensembl.org/overlap/region/human/{ch}:{pos}-{pos}?feature=variation;content-type=application/json'
  j=get(url,f'{ch}.{pos}.overlap.json');ids=sorted(set(v['id'] for v in j if v.get('assembly_name')=='GRCh38' and v.get('start')==pos and v.get('end')==pos and ref in v.get('alleles',[]) and alt in v.get('alleles',[])))
  r['exact_public_ids']=ids;r['exact_forward_mappings']=[];r['global_alt_frequencies']=[];r['variation_metadata_significance_not_allele_specific']=[]
  for vid in ids:
   u=f'https://rest.ensembl.org/variation/human/{vid}?pops=1;content-type=application/json';v=get(u,vid+'.variation.json')
   maps=[m for m in v.get('mappings',[]) if m.get('assembly_name')=='GRCh38' and m.get('seq_region_name')==ch and m.get('start')==pos and m.get('end')==pos and m.get('strand')==1 and ref in m.get('allele_string','').split('/') and alt in m.get('allele_string','').split('/')]
   if not maps:raise ValueError('Exact forward mapping unavailable '+vid)
   r['exact_forward_mappings']+=maps;r['variation_metadata_significance_not_allele_specific']+=v.get('clinical_significance',[])
   for f in v.get('populations',[]):
    if f.get('allele')==alt and f.get('frequency') is not None and f.get('population') in {'gnomADe:ALL','gnomADg:ALL','1000GENOMES:phase_3:ALL'}:r['global_alt_frequencies'].append({'id':vid,**f})
  r['status']='exact_catalog_match' if ids else 'no_exact_match_returned'
 except Exception as e:r.update(status='lookup_error',error=type(e).__name__+': '+str(e))
 print(c['gene'],c['protein'],r['status'],r.get('exact_public_ids'),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:pub=list(pool.map(one,[c for c in C if c['consequence']=='missense']))
o={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'17 missense variants queried by generic public genomic locus/allele only, without patient identifier or raw sequencing data. Exact GRCh38 coordinate, ref, alt, and forward orientation verified before alt-frequency use. Catalog metadata is not clinical adjudication or constitutional status.','VCF_exact_comparison':all_exact,'missense_public_context':pub}
(B/'independent-public-audit.json').write_text(json.dumps(o,indent=2)+'\n')

"""Public normal catalog lookup only; not germline or somatic classification."""
import pathlib,csv,json,urllib.request,urllib.error,time,collections,datetime
B=pathlib.Path(__file__).resolve().parent;ROOT=B.parents[2]
C=[{'gene':r['gene'],'protein':'p.'+r['protein'],'hg38_chrom':r['chrom'],'position_1based':str(r['pos1']),'ref':r['ref'],'alt':r['alt'],'DNA_allele_fraction':(r['Caris_VCF_exact_matches'][0]['VAF'] if r['Caris_VCF_exact_matches'] else str(r['screen_fraction'])),'DNA_allele_fraction_source':('Caris VCF' if r['Caris_VCF_exact_matches'] else 'research read-count screen')} for r in json.loads((B/'screen-results.json').read_text())['candidates'] if r['consequence']=='missense'];assert len(C)==12
CACHE=B/'public-population';CACHE.mkdir(exist_ok=True)
def get(url,file):
 p=CACHE/file
 if p.exists():return json.loads(p.read_text())['result']
 for n in range(4):
  try:
   req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'Caris-research-public-reference-annotation/1.0'})
   with urllib.request.urlopen(req,timeout=45) as r:x=json.load(r)
   p.write_text(json.dumps({'source_url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'result':x},indent=2)+'\n');return x
  except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError):
   if n==3:raise
   time.sleep(1+n)
rows=[];detail=[]
for c in C:
 gene=c['gene'];pos=int(c['position_1based']);ch=c['hg38_chrom'].replace('chr','');ref=c['ref'];alt=c['alt'];key=f'{gene}.{pos}'
 r={'gene':gene,'protein':c['protein'],'chrom':c['hg38_chrom'],'position_1based':pos,'ref':ref,'alt':alt,'DNA_allele_fraction':c['DNA_allele_fraction'],'DNA_allele_fraction_source':c['DNA_allele_fraction_source'],'status':'queried'}
 try:
  url=f'https://rest.ensembl.org/overlap/region/human/{ch}:{pos}-{pos}?feature=variation;content-type=application/json'
  overlaps=get(url,key+'.overlap.json')
  exact=[v for v in overlaps if v.get('assembly_name')=='GRCh38' and v.get('start')==pos and v.get('end')==pos and ref in v.get('alleles',[]) and alt in v.get('alleles',[])]
  ids=sorted(set(v['id'] for v in exact));r['public_exact_variant_ids']=';'.join(ids)
  infos=[];freqs=[];sigs=[];synonyms=[];mappings=[]
  for vid in ids:
   vurl=f'https://rest.ensembl.org/variation/human/{vid}?pops=1;content-type=application/json'
   j=get(vurl,vid+'.variation.json');infos.append({'id':vid,'url':vurl,'result':j})
   same=[m for m in j.get('mappings',[]) if m.get('assembly_name')=='GRCh38' and m.get('seq_region_name')==ch and m.get('start')==pos and m.get('end')==pos and m.get('strand')==1 and ref in m.get('allele_string','').split('/') and alt in m.get('allele_string','').split('/')]
   if not same:raise RuntimeError('No exact forward-strand GRCh38 mapping for '+vid)
   mappings+=same;sigs+=j.get('clinical_significance',[]);synonyms+=j.get('synonyms',[])
   for q in j.get('populations',[]):
    if q.get('allele')==alt and q.get('frequency') is not None:
     freqs.append({'id':vid,'population':q['population'],'allele':alt,'frequency':q['frequency']})
  allfreq=[f for f in freqs if f['population'] in {'gnomADe:ALL','gnomADg:ALL','1000GENOMES:phase_3:ALL','ESP6500:ALL','TOPMED'}]
  def freq(pop):
   a=[float(f['frequency']) for f in freqs if f['population']==pop];return max(a) if a else ''
  for pop in ['gnomADe:ALL','gnomADg:ALL','1000GENOMES:phase_3:ALL']:
   r[pop+'_observed_alt_frequency']=freq(pop)
  gnom=[f for f in freqs if f['population'].startswith(('gnomADe:','gnomADg:'))];r['gnomAD_any_population_max_alt_frequency']=max([float(f['frequency']) for f in gnom],default='')
  r['Ensembl_clinical_significance_not_allele_specific']=';'.join(sorted(set(sigs)))
  r['clinical_metadata_warning']='Clinical significance is public variation-record metadata, potentially multi-allelic/condition-dependent; not clinical adjudication of this specific allele.'
  r['public_synonym_ids']=';'.join(sorted(set(s for s in synonyms if s.startswith(('VCV','RCV','COSM','COSV')))))
  major=[r[pop+'_observed_alt_frequency'] for pop in ['gnomADe:ALL','gnomADg:ALL','1000GENOMES:phase_3:ALL'] if r[pop+'_observed_alt_frequency']!='']
  if major and max(major)>=.01:r['neoantigen_triage']='Deprioritize as tumor-specific: observed allele common (>=1% in a global reference population); germline test still establishes patient status.'
  elif major and max(major)>0:r['neoantigen_triage']='Rare population allele present; not proven somatic. Matched normal required before tumor-specific targeting.'
  elif ids:r['neoantigen_triage']='Public exact variant exists; population frequency unavailable or zero in returned datasets. Somatic status unresolved.'
  else:r['neoantigen_triage']='No exact match returned by limited Ensembl catalog lookup; not proof of novelty or somatic status.'
  detail.append({'candidate':r,'overlap_url':url,'exact_overlaps':exact,'exact_mappings':mappings,'alt_population_frequencies':freqs,'global_alt_population_frequencies':allfreq,'variation_records_cached':[j['id'] for j in infos]})
 except Exception as e:r.update(status='lookup_error',error=type(e).__name__+': '+str(e));detail.append({'candidate':r})
 rows.append(r);print(gene,c['protein'],r.get('public_exact_variant_ids'),r.get('gnomADe:ALL_observed_alt_frequency'),r.get('neoantigen_triage'),flush=True)
 (B/'extended-protein-candidate-public-annotation.checkpoint.json').write_text(json.dumps(detail,indent=2)+'\n')
keys=list(dict.fromkeys(k for r in rows for k in r))
with (B/'extended-protein-candidate-public-annotation.tsv').open('w') as f:
 w=csv.DictWriter(f,keys,delimiter='\t');w.writeheader();w.writerows(rows)
(B/'extended-protein-candidate-public-annotation.json').write_text(json.dumps({'scope':'12 immune-presentation panel protein-changing candidates; public Ensembl/dbSNP population catalog reference queries only. No patient IDs, read sequences, raw sequencing data or patient files transmitted.','limitations':['Population presence does not prove this patient constitutional status; absent catalog allele does not establish somatic status.','Frequency is allele-specific only after exact GRCh38 forward-strand coordinate/ref/alt mapping verification.','Ensembl clinical-significance field may combine multiple alternate alleles and conditions; not used to assign pathogenicity.','Public datasets may change and coverage/ancestry representation is uneven; no frequencies treated as patient probabilities.'],'results':detail},indent=2)+'\n')

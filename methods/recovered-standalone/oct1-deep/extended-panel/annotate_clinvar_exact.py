"""Exact GRCh38 REF/ALT ClinVar lookup from public catalog VCV links."""
import pathlib,json,csv,urllib.request,urllib.error,time,xml.etree.ElementTree as ET,datetime
B=pathlib.Path(__file__).resolve().parent;CACHE=B/'public-population'
x=json.loads((B/'extended-protein-candidate-public-annotation.json').read_text());out=[]
for item in x['results']:
 c=item['candidate'];row={k:c.get(k,'') for k in ['gene','protein','chrom','position_1based','ref','alt','public_exact_variant_ids']};matches=[];errors=[]
 vcvs=[s for s in c.get('public_synonym_ids','').split(';') if s.startswith('VCV')]
 for vcv in vcvs:
  file=CACHE/f'{vcv}.clinvar.xml';url=f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&id={vcv}'
  try:
   if not file.exists():
    time.sleep(.4)
    for attempt in range(3):
     try:
      with urllib.request.urlopen(url,timeout=40) as r:raw=r.read()
      ET.fromstring(raw);file.write_bytes(raw);break
     except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError):
      if attempt==2:raise
      time.sleep(2+attempt)
   xml=ET.parse(file).getroot()
   for va in xml.findall('.//VariationArchive'):
    locs=[s.attrib for s in va.findall('.//SequenceLocation') if s.get('Assembly')=='GRCh38' and s.get('Chr')==str(c['chrom']).replace('chr','') and s.get('positionVCF')==str(c['position_1based']) and s.get('referenceAlleleVCF')==c['ref'] and s.get('alternateAlleleVCF')==c['alt']]
    if not locs:continue
    cr=va.find('ClassifiedRecord');m={'accession':va.get('Accession'),'version':va.get('Version'),'name':va.get('VariationName'),'variation_id':va.get('VariationID'),'last_updated':va.get('DateLastUpdated'),'source_url':url,'browser_url':f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{va.get("VariationID")}/','exact_GRCh38_locations':locs}
    for label in ['GermlineClassification','SomaticClinicalImpact','OncogenicityClassification']:
     el=cr.find('Classifications/'+label) if cr is not None else None
     m[label]=None if el is None else {'description':el.findtext('Description'),'review_status':el.findtext('ReviewStatus'),'attributes':el.attrib,'other_text':' '.join(z.strip() for z in el.itertext() if z.strip())}
    matches.append(m)
  except Exception as e:errors.append({'vcv':vcv,'error':type(e).__name__+': '+str(e)})
 row['exact_ClinVar_accessions']=';'.join(m['accession']+'.'+m['version'] for m in matches)
 row['ClinVar_germline_aggregate']=';'.join(m['GermlineClassification']['description'] or '' for m in matches if m['GermlineClassification'])
 row['ClinVar_germline_review_status']=';'.join(m['GermlineClassification']['review_status'] or '' for m in matches if m['GermlineClassification'])
 row['ClinVar_oncogenicity_aggregate']=';'.join(m['OncogenicityClassification']['description'] or '' for m in matches if m['OncogenicityClassification'])
 row['ClinVar_somatic_clinical_impact']=';'.join(m['SomaticClinicalImpact']['description'] or '' for m in matches if m['SomaticClinicalImpact'])
 row['ClinVar_links']=';'.join(m['browser_url'] for m in matches)
 row['status']='exact matched' if matches else 'lookup error' if errors else 'no exact linked ClinVar VCV returned'
 row['interpretation_limit']='Germline pathogenicity classifications do not prove patient germline status, somatic oncogenicity, immunogenicity or therapeutic actionability. No-entry is not benign or pathogenic.'
 out.append({'row':row,'matches':matches,'errors':errors});print(c['gene'],c['protein'],row['exact_ClinVar_accessions'],row['ClinVar_germline_aggregate'],flush=True)
 (B/'extended-protein-candidate-clinvar.checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
with (B/'extended-protein-candidate-clinvar.tsv').open('w') as f:
 w=csv.DictWriter(f,list(out[0]['row']),delimiter='\t');w.writeheader();w.writerows(x['row'] for x in out)
(B/'extended-protein-candidate-clinvar.json').write_text(json.dumps({'method':'Retrieve public VCVs linked in exact-allele Ensembl catalog records; retain only exact forward-genome GRCh38 chromosome/positionVCF/referenceAlleleVCF/alternateAlleleVCF matches from NCBI ClinVar XML. Public identifiers only; no patient IDs or sequences uploaded.','documentation':'https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/','retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':out},indent=2)+'\n')

"""Public catalog lookup by region, with exact haplotype matching locally.
Only generic public genomic-region/dbSNP requests; no patient identity, files,
sequences, composite variant list or medical history transmitted.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,urllib.request,urllib.error,time,datetime,csv
from Bio.Seq import Seq
OUT=Path(__file__).resolve().parent; CACHE=OUT/'public-catalog';CACHE.mkdir(exist_ok=True)
CAND=json.loads((OUT/'expressed-frameshift-shortlist.json').read_text())['candidates']

def get(url,name):
    p=CACHE/name
    if p.exists():return json.loads(p.read_text())['result']
    for a in range(4):
        try:
            req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'Public-reference-research/1.0'})
            with urllib.request.urlopen(req,timeout=50) as r:data=json.load(r)
            p.write_text(json.dumps({'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'result':data},indent=2)+'\n');return data
        except Exception:
            if a==3:raise
            time.sleep(a+1)

def matching_alts(v,start,end,alleles,window):
    """Ensembl 1-based inclusive intervals, including insertion start>end."""
    ws=window['start'];seq=window['dna'].upper();s=start-1;e=end
    if s<ws or e>window['end'] or s>e:return []
    # An insertion has start=end+1 -> s==e, which is correct here.
    original=seq[s-ws:e-ws]
    if original not in [a.replace('-','') for a in alleles]:return []
    p=v['pos1']-1;target=seq[:p-ws]+v['alt']+seq[p-ws+len(v['ref']):]
    return [a for a in alleles if seq[:s-ws]+a.replace('-','')+seq[e-ws:]==target]

def annotate(v):
    chrom=v['chrom'].removeprefix('chr');p=v['pos1']-1
    ws=max(0,p-80);we=p+len(v['ref'])+80
    window=json.loads((OUT/'reference-windows'/f'{v["chrom"]}-{ws}-{we}.json').read_text())
    start=max(1,v['pos1']-40);end=v['pos1']+len(v['ref'])+40
    url=f'https://rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}?feature=variation;content-type=application/json'
    result={'gene':';'.join(v['genes']),'chrom':v['chrom'],'pos1':v['pos1'],'ref':v['ref'],'alt':v['alt'],'DNA_fraction':v['DNA']['alt_fraction_ref_alt_only'],'RNA_alt_names':v['RNA']['paired_names'].get('alt',0),'matches':[],'lookup_status':'success'}
    try:
        overlaps=get(url,f'{chrom}-{start}-{end}.overlap.json')
        ids=set()
        for q in overlaps:
            if q.get('assembly_name')!='GRCh38':continue
            alleles=q.get('alleles',[])
            if q.get('strand',1)==-1:alleles=[str(Seq(a).reverse_complement()) if a!='-' else a for a in alleles]
            if matching_alts(v,q['start'],q['end'],alleles,window):ids.add(q['id'])
        for ident in sorted(ids):
            data=get(f'https://rest.ensembl.org/variation/human/{ident}?pops=1;content-type=application/json',f'{ident}.variation.json')
            mapped=[];popalts=set()
            for m in data.get('mappings',[]):
                if m.get('assembly_name')!='GRCh38' or m.get('seq_region_name')!=chrom or m.get('strand')!=1:continue
                same=matching_alts(v,m['start'],m['end'],m.get('allele_string','').split('/'),window)
                if same:mapped.append(m);popalts.update(same)
            if not mapped:continue
            freq=[q for q in data.get('populations',[]) if q.get('allele') in popalts and q.get('frequency') is not None]
            globalfreq=[q for q in freq if q['population'] in ['gnomADe:ALL','gnomADg:ALL','1000GENOMES:phase_3:ALL','TOPMED','ESP6500:ALL']]
            result['matches'].append({'id':ident,'matched_mappings':mapped,'matched_population_alleles':sorted(popalts),'global_frequencies':globalfreq,'synonyms':data.get('synonyms',[]),'clinical_metadata_not_used_for_allele_classification':data.get('clinical_significance',[])})
        afs=[float(f['frequency']) for m in result['matches'] for f in m['global_frequencies']]
        result['largest_returned_global_alt_frequency']=max(afs) if afs else None
        result['interpretation']=('Common catalog allele (>=1% global frequency), not a tumor-specific vaccine target without contrary matched-normal evidence.' if afs and max(afs)>=.01 else 'Catalog allele present; inherited versus somatic status in this patient unresolved.' if result['matches'] else 'No exact public haplotype match returned; not evidence of novelty or somatic status.')
    except Exception as e:result['lookup_status']='error';result['error']=repr(e)
    print(result['gene'],[m['id'] for m in result['matches']],result.get('largest_returned_global_alt_frequency'),flush=True)
    return result

if __name__=='__main__':
    rows=list(ThreadPoolExecutor(max_workers=4).map(annotate,CAND))
    (OUT/'expressed-indel-public-context.json').write_text(json.dumps({'method':'Public Ensembl/dbSNP overlap and variation lookups; exact local GRCh38 reference/alternate haplotype equivalence, then forward-strand mapping validation before frequency attribution. Rare/absent catalog alleles remain unclassified. Returned frequency datasets are version-specific and may have uneven indel coverage. No clinical significance assigned from record-wide metadata.','results':rows},indent=2)+'\n')
    fields=['gene','chrom','pos1','ref','alt','DNA_fraction','RNA_alt_names','public_ids','largest_returned_global_alt_frequency','lookup_status','interpretation']
    with (OUT/'expressed-indel-public-context.tsv').open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader()
        for r in rows:w.writerow({**{k:r.get(k,'') for k in fields},'public_ids':';'.join(m['id'] for m in r['matches'])})

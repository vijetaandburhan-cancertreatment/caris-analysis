"""Only fixed public references; no patient data are sent."""
import concurrent.futures,datetime,hashlib,json,pathlib,time,urllib.request
OUT=pathlib.Path(__file__).resolve().parent/'public'
BASE='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/'
URLS={
 '1000G.GRCh38.20190312.sites.vcf.gz':BASE+'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz',
 '1000G.GRCh38.20190312.sites.vcf.gz.tbi':BASE+'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz.tbi',
 '1000G.README.txt':BASE+'20190312_biallelic_SNV_and_INDEL_README.txt',
 '1000G.MANIFEST.txt':BASE+'20190312_biallelic_SNV_and_INDEL_MANIFEST.txt',
 'hg38.2bit':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit'
}
def get(item):
 name,url=item;p=OUT/name;started=time.time();h=hashlib.sha256();n=0
 if p.exists():
  with p.open('rb') as f:
   for z in iter(lambda:f.read(1024*1024),b''):h.update(z);n+=len(z)
  return {'file':str(p),'url':url,'bytes':n,'sha256':h.hexdigest(),'reused':True}
 with urllib.request.urlopen(url,timeout=120) as r,(OUT/(name+'.partial')).open('wb') as f:
  last=time.time()
  for z in iter(lambda:r.read(1024*1024),b''):
   f.write(z);h.update(z);n+=len(z)
   if time.time()-last>20:print(name,n,round(time.time()-started,1),flush=True);last=time.time()
 (OUT/(name+'.partial')).rename(p)
 res={'file':str(p),'url':url,'bytes':n,'sha256':h.hexdigest(),'elapsed_seconds':time.time()-started,'downloaded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};print(name,'complete',n,flush=True);return res
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:results=list(ex.map(get,URLS.items()))
(OUT/'manifest.json').write_text(json.dumps(results,indent=2)+'\n')

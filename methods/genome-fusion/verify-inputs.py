import pathlib,json,hashlib,datetime,gzip,itertools
P=pathlib.Path(__file__).resolve().parent;old=json.loads((P.parents[1]/'oct1-analysis/resident-copies-verification.json').read_text());prior={x['name']:x for x in old['files']};raw=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');out={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':[],'originals_modified':False}
paths=[]
for mate in [1,2]:
 name=f'RNA_TN26-279853_S25.R{mate}.fastq.gz';f=raw/name;paths.append(f);h=hashlib.sha256();s=f.stat()
 with f.open('rb') as stream:
  for chunk in iter(lambda:stream.read(8*1024**2),b''):h.update(chunk)
 now=f.stat();entry={'path':str(f),'size_bytes':s.st_size,'mtime_ns':s.st_mtime_ns,'sha256':h.hexdigest(),'matches_previously_source_CRC_verified_SHA256':h.hexdigest()==prior[name]['calculated']['sha256'],'unchanged_during_read':(s.st_size,s.st_mtime_ns)==(now.st_size,now.st_mtime_ns)};assert entry['matches_previously_source_CRC_verified_SHA256'] and entry['unchanged_during_read'];out['files'].append(entry)
errors=[];n=0;maxlen=0
with gzip.open(paths[0],'rb') as a,gzip.open(paths[1],'rb') as b:
 for i in range(10000):
  x=[a.readline().rstrip(b'\r\n') for _ in range(4)];y=[b.readline().rstrip(b'\r\n') for _ in range(4)]
  if not x[0] or not y[0]:errors.append('unexpected early EOF');break
  names=[z[0].split()[0].removesuffix(b'/1').removesuffix(b'/2') for z in [x,y]]
  if names[0]!=names[1]:errors.append(f'pair mismatch{i}')
  for z in [x,y]:
   if not z[0].startswith(b'@') or not z[2].startswith(b'+') or len(z[1])!=len(z[3]):errors.append(f'format error{i}')
   maxlen=max(maxlen,len(z[1]))
  n+=1
out['paired_preflight']={'first_pairs_checked':n,'errors':errors,'maximum_read_length_in_sample':maxlen,'scope':'First10000 pairs only; exact whole-file SHA256 matches prior source-CRC-verified files. Whole-file processed count expected23209264 pairs from prior full kallisto run; STAR final count will be checked independently.'};assert not errors;out['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(P/'input-reverification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

import pathlib,json,gzip,hashlib,collections,time,resource,shutil,pysam,ahocorasick_rs
O=pathlib.Path(__file__).parent;B=pathlib.Path.home()/'.local/share/codex/caris-analysis';rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def guard():
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*2**30:raise RuntimeError('RSS limit4GiB')
 if shutil.disk_usage(O).free<3*2**30:raise RuntimeError('Free-space guard3GiB')
guard();t0=time.time();q=json.loads((O/'patient-queries.json').read_text());controls=json.loads((O/'controls.json').read_text());q+=controls
genome=B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa';transcripts=B/'oct4-hla-allele-support/gencode.v37.transcripts.fa.gz'
oriented=[]
for x in q:
 for orientation,seq in[('+',x['sequence']),('-',rc(x['sequence']))]:oriented.append(dict(x,oriented_id=x['query_id']+orientation,orientation=orientation,oriented_sequence=seq))
seedroles=collections.defaultdict(list)
for x in oriented:
 s=x['oriented_sequence'];starts=list(range(0,len(s)-18,19));starts=sorted(set(starts+[len(s)-19]));x['seed_starts']=starts
 for st in starts:
  seed=s[st:st+19]
  if set(seed)<=set('ACGT'):seedroles[seed].append((x['oriented_id'],st))
patterns=sorted(seedroles);auto=ahocorasick_rs.AhoCorasick(patterns);hits={'genome':collections.defaultdict(list),'transcript':collections.defaultdict(list)};counts={'genome':collections.Counter(),'transcript':collections.Counter()};caps={'genome':set(),'transcript':set()};cap=1000
# Control seed matcher against independent plain substring search, including overlapping repeat hits.
for s in [controls[0]['sequence'],controls[-1]['sequence'],controls[2]['sequence']]:
 got=sorted(auto.find_matches_as_indexes(s,overlapping=True));want=sorted((i,j,j+len(p))for i,p in enumerate(patterns)for j in range(len(s)-len(p)+1)if s[j:j+len(p)]==p);assert got==want

def scan(kind,refid,s,offset=0,tail_length=0):
 for i,a,b in auto.find_matches_as_indexes(s,overlapping=True):
  if b<=tail_length:continue
  counts[kind][i]+=1
  if len(hits[kind][i])<cap:hits[kind][i].append((refid,offset+a))
  else:caps[kind].add(i)
# streaming chunk genome;18nt overlap preserves every19mer exactly once.
nbase=ncontig=0;contig=None;buf=[];n=0;tail='';consumed=0
with genome.open()as f:
 for line in f:
  if line.startswith('>'):
   if buf:
    s=''.join(buf);scan('genome',contig,tail+s,consumed-len(tail),len(tail));buf=[];n=0
   contig=line[1:].split()[0];ncontig+=1;tail='';consumed=0
  else:
   s=line.strip().upper();buf.append(s);n+=len(s);nbase+=len(s)
   if n>=1048576:
    s=''.join(buf);scan('genome',contig,tail+s,consumed-len(tail),len(tail));consumed+=len(s);tail=(tail+s)[-18:];buf=[];n=0;guard()
 if buf:s=''.join(buf);scan('genome',contig,tail+s,consumed-len(tail),len(tail))
print('genome complete',nbase,time.time()-t0,flush=True)
ntx=txbase=0
with pysam.FastxFile(str(transcripts))as f:
 for tx in f:
  scan('transcript',tx.name.split('|')[0],tx.sequence.upper());ntx+=1;txbase+=len(tx.sequence)
  if ntx%10000==0:guard()
# Nominate origins with20nt slop, before scoring. Merge only origin positions within5nt to preserve distinct placements.
byquery={x['oriented_id']:{'genome':collections.defaultdict(set),'transcript':collections.defaultdict(set)}for x in oriented}
for kind in hits:
 for i,ls in hits[kind].items():
  for oid,qpos in seedroles[patterns[i]]:
   for refid,pos in ls:byquery[oid][kind][refid].add(pos-qpos)
outputs=[];windowcaps=[]
for x in oriented:
 v={'query':x,'references':{}}
 for kind,rrs in byquery[x['oriented_id']].items():
  ww=[]
  for refid,starts in rrs.items():
   ss=sorted(starts);groups=[]
   for st in ss:
    if groups and st-groups[-1][-1]<=5:groups[-1].append(st)
    else:groups.append([st])
   for g in groups:ww.append({'reference_id':refid,'start0':max(0,g[0]-20),'end0':max(0,g[-1]+len(x['oriented_sequence'])+20),'seed_origin_min':g[0],'seed_origin_max':g[-1],'seed_origin_count':len(g)})
  if len(ww)>500:windowcaps.append({'query':x['oriented_id'],'reference':kind,'nominated':len(ww)});ww=ww[:500]
  v['references'][kind]=ww
 outputs.append(v)
report={'queries':outputs,'controls':controls,'seed_patterns':patterns,'seed_roles':{str(i):seedroles[s]for i,s in enumerate(patterns)},'seed_hit_counts':{k:dict(v)for k,v in counts.items()},'seed_caps':{k:sorted(v)for k,v in caps.items()},'window_caps':windowcaps,'genome_bases':nbase,'genome_contigs':ncontig,'transcript_bases':txbase,'transcripts':ntx,'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'exact_seed_public_regression_pass':True}
(O/'nominations.json').write_text(json.dumps(report,indent=2)+'\n');(O/'queries-with-controls.json').write_text(json.dumps(q,indent=2)+'\n')
print(json.dumps({'controls':len(controls),'oriented_queries':len(outputs),'seed_patterns':len(patterns),'windows':{k:sum(len(x['references'][k])for x in outputs)for k in hits},'capped_seeds':{k:len(v)for k,v in caps.items()},'window_caps':windowcaps,'elapsed_seconds':report['elapsed_seconds'],'peak_RSS_bytes':report['peak_RSS_bytes']}))

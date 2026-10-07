import pathlib,json,gzip,hashlib,pysam,itertools,collections,shutil
B=pathlib.Path.home()/'.local/share/codex/caris-analysis';A=B/'oct4-fusion-read-audit';O=pathlib.Path(__file__).parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
d=json.loads((A/'patient/junction-read-audit.json').read_text());r=next(x for x in d['rows']if x['candidate_id']=='cb8931f687fc');names=set(r['read_names']);marker=next(x['sequence']for x in r['junction_patterns']if x['arm_bases']==25);out=[]
with pysam.FastxFile(str(A/'patient/audit-selected.R1.fastq.gz'))as a,pysam.FastxFile(str(A/'patient/audit-selected.R2.fastq.gz'))as b:
 for r1,r2 in itertools.zip_longest(a,b):
  assert r1 and r2 and r1.name==r2.name;hits=[]
  for mate,read in enumerate((r1,r2),1):
   for ori,pat in [('+',marker),('-',rc(marker))]:
    for k in range(len(read.sequence)-49):
     if read.sequence[k:k+50]==pat:
      q=min(ord(x)-33for x in read.quality[k:k+50]);hits.append({'mate':mate,'orientation':ori,'start0':k,'Qmin':q})
  isexact=any(x['Qmin']>=20for x in hits)
  if r1.name in names or isexact:
   out.append({'name':r1.name,'name_sha256':hashlib.sha256(r1.name.encode()).hexdigest(),'caller_named':r1.name in names,'Q20_50nt_support':isexact,'hits':hits,'R1_sequence':r1.sequence,'R1_quality':r1.quality,'R2_sequence':r2.sequence,'R2_quality':r2.quality,'orientation_normalized_pair':min(r1.sequence+'|'+r2.sequence,rc(r2.sequence)+'|'+rc(r1.sequence))})
assert sum(x['Q20_50nt_support']for x in out)==11
for n,x in enumerate(out,1):x['pair_id']=f'P{n:02d}'
family={s:i+1 for i,s in enumerate(sorted({x['orientation_normalized_pair']for x in out}))}
for x in out:x['sequence_family']=f'F{family[x.pop("orientation_normalized_pair")]:02d}'
querygroups=collections.defaultdict(list)
for x in out:
 for m in(1,2):querygroups[x[f'R{m}_sequence']].append({'pair_id':x['pair_id'],'mate':m,'quality':x[f'R{m}_quality'],'caller_named':x['caller_named'],'Q20_50nt_pair_support':x['Q20_50nt_support'],'family':x['sequence_family']})
qs=[{'query_id':f'Q{i:02d}','sequence':s,'members':ms,'kind':'patient'}for i,(s,ms)in enumerate(sorted(querygroups.items()),1)]
# Load selected parent transcripts and reference-derived controls, no patient-data fetch.
wanted={'ENST00000311127.9','ENST00000397752.8'};ts=[]
with pysam.FastxFile(str(B/'oct4-hla-allele-support/gencode.v37.transcripts.fa.gz'))as f:
 for x in f:
  fields=x.name.split('|')
  if fields[0]in wanted or (len(fields)>5 and fields[5]=='CCND1'):ts.append({'id':fields[0],'name':x.name,'sequence':x.sequence.upper()})
(O/'selected-parent-transcripts.json').write_text(json.dumps(ts,indent=2)+'\n')
(O/'pairs.json').write_text(json.dumps(out,indent=2)+'\n');(O/'patient-queries.json').write_text(json.dumps(qs,indent=2)+'\n')
r.pop('read_names',None);r.pop('read_identifiers',None);(O/'source-candidate.json').write_text(json.dumps(r,indent=2)+'\n')
prov={'prior_raw_audit_sha256':sha(A/'patient/junction-read-audit.json'),'selected_R1_sha256':sha(A/'patient/audit-selected.R1.fastq.gz'),'selected_R2_sha256':sha(A/'patient/audit-selected.R2.fastq.gz'),'pair_count':len(out),'caller_named_pairs':sum(x['caller_named']for x in out),'Q20_exact50_pairs':sum(x['Q20_50nt_support']for x in out),'whole_pair_sequence_families_all':len(family),'whole_pair_sequence_families_Q20':len({x['sequence_family']for x in out if x['Q20_50nt_support']}),'distinct_mate_sequences':len(qs),'originals_modified':False}
(O/'preparation-provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
with gzip.open(A/'patient/selected-Chimeric.out.junction.gz','rt')as f,(O/'selected-STAR-Chimeric.out.junction').open('w')as g:
 for line in f:
  z=line.split('\t')
  if line.startswith('#')or(len(z)>9 and z[9]in{x['name']for x in out}):g.write(line)
print(json.dumps(prov));print('parents',[(x['id'],len(x['sequence']))for x in ts]);print('query lengths',collections.Counter(len(x['sequence'])for x in qs))

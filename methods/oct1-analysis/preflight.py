from pathlib import Path
import gzip,struct,json,collections,hashlib,datetime
p=Path('outputs/caris-raw-data/TN26-279853');out={}
eof=bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
for name in ['DNA_TN26-279853.bam','RNA_TN26-279853.bam']:
 f=p/name
 with f.open('rb') as h:h.seek(-28,2);has_eof=h.read()==eof
 with gzip.open(f,'rb') as h:
  assert h.read(4)==b'BAM\1'
  l=struct.unpack('<i',h.read(4))[0];header=h.read(l).decode();n=struct.unpack('<i',h.read(4))[0]
  refs=[]
  for _ in range(n):
   ln=struct.unpack('<i',h.read(4))[0];r=h.read(ln).rstrip(b'\0').decode();length=struct.unpack('<i',h.read(4))[0];refs.append({'name':r,'length':length})
 rows=collections.defaultdict(list)
 for line in header.splitlines():
  pieces=line.split('\t');rows[pieces[0]].append(dict(x.split(':',1) for x in pieces[1:] if ':' in x))
 Path('work/oct1-analysis/'+name+'.header.sam').write_text(header)
 out[name]={'size_bytes':f.stat().st_size,'bam_bgzf_eof_present':has_eof,'HD':rows['@HD'],'read_groups':rows['@RG'],'programs':[{k:v for k,v in z.items() if k!='CL'} for z in rows['@PG']],'program_command_summary':[{'ID':z.get('ID'),'CL':z.get('CL','')[:1500]} for z in rows['@PG']],'reference_count':n,'references_first_30':refs[:30],'reference_header_first_3':rows['@SQ'][:3],'comments':rows['@CO']}
Path('work/oct1-analysis/bam-header-preflight.json').write_text(json.dumps(out,indent=2))
for prefix in ['DNA']:
 fs=[gzip.open(p/(prefix+'_TN26-279853_S25.R'+str(i)+'.fastq.gz'),'rt') for i in (1,2)]
 counts=[];mismatches=0;lengths=[collections.Counter(),collections.Counter()];q=[0,0];q30=[0,0];bases=[0,0];head=[]
 for i in range(10000):
  rs=[]
  for j,f in enumerate(fs):
   lines=[f.readline().rstrip('\n') for _ in range(4)]
   if not lines[0]:break
   assert lines[0].startswith('@') and lines[2].startswith('+') and len(lines[1])==len(lines[3])
   rs.append(lines);lengths[j][len(lines[1])]+=1;bases[j]+=len(lines[1]);q[j]+=sum(ord(c)-33 for c in lines[3]);q30[j]+=sum(ord(c)-33>=30 for c in lines[3])
  if len(rs)!=2:break
  a,b=[z[0].split()[0].removesuffix('/1').removesuffix('/2') for z in rs]
  mismatches+=a!=b
  if i==0:head=[z[0] for z in rs]
 for f in fs:f.close()
 out[prefix+'_FASTQ_first_10000_pairs']={'pairs_sampled':i+1,'mate_name_mismatches':mismatches,'read_length_counts':[dict(z) for z in lengths],'mean_Q':[q[j]/bases[j] for j in (0,1)],'fraction_Q30':[q30[j]/bases[j] for j in (0,1)],'first_read_headers':head,'sampling_limit':'First10,000pairs only; notrepresentativefull-fileQC'}
Path('work/oct1-analysis/preflight.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))

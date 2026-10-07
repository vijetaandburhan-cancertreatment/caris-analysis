#!/usr/bin/env python3
import collections,gzip,json,pathlib
base=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-fusion-read-audit');a=json.load(open(base/'patient/junction-read-audit.json'));rows=[r for r in a['rows'] if r['source_set']=='accepted']; byname=collections.defaultdict(list)
def ep(bp,d):c,p=bp.split(':');return(c,d,int(p))
def order(*ep):return sorted(ep)
for r in rows:
 r['_ep']=order(*(ep(r[f'breakpoint{i}'],r[f'direction{i}']) for i in [1,2]))
 for n in r['read_names']:byname[n].append(r['audit_key'])
rmap={r['audit_key']:r for r in rows};stats={r['audit_key']:collections.defaultdict(set) for r in rows};details={r['audit_key']:collections.defaultdict(list) for r in rows}
with gzip.open(base/'patient/selected-Chimeric.out.junction.gz','rt') as f:
 for line in f:
  z=line.rstrip().split('\t')
  if len(z)<14 or not z[1].isdigit() or z[9] not in byname:continue
  endpoints=order((z[0],'downstream' if z[2]=='+' else 'upstream',int(z[1])-1 if z[2]=='+' else int(z[1])+1),(z[3],'upstream' if z[5]=='+' else 'downstream',int(z[4])+1 if z[5]=='+' else int(z[4])-1))
  for cid in byname[z[9]]:
   r=rmap[cid];s=stats[cid];s['names_with_any_chimeric_record'].add(z[9]);kind='other_breakends'
   if all(x[:2]==y[:2] for x,y in zip(endpoints,r['_ep'])):
    delta=[x[2]-y[2] for x,y in zip(endpoints,r['_ep'])]
    if all(d==0 for d in delta):kind='exact_reported_breakends'
    elif all(abs(d)<=25 for d in delta):kind='within25bp_breakends'
   s[kind].add(z[9]);details[cid][z[9]].append({'kind':kind,'endpoints':endpoints,'num_chim_aln':int(z[14]) if len(z)>14 else None,'scores':z[15:19]})
   if len(z)>14 and int(z[14])>1:s['names_with_multiple_STAR_chimeric_alignments'].add(z[9])
output=[]
for r in rows:
 cid=r['audit_key'];d={'candidate_id':cid,'gene1':r['gene1'],'gene2':r['gene2'],'caller_named_pairs':len(r['read_names'])}
 for k in ['names_with_any_chimeric_record','exact_reported_breakends','within25bp_breakends','other_breakends','names_with_multiple_STAR_chimeric_alignments']:d[k]=len(stats[cid][k])
 d['names_no_chimeric_record']=len(set(r['read_names'])-stats[cid]['names_with_any_chimeric_record']);output.append(d)
(base/'patient-evidence/STAR-name-diagnostics.json').write_text(json.dumps({'rows':output,'name_level_records':details,'limits':['Endpoint proximity is not Arriba filtered read membership. Distinct alternatives reflect STAR records and do not establish biological origin.','STAR merged paired-end CIGARs can represent a fragment instead of a single sequenced read.']},indent=2)+'\n')
for r in output:print(r)

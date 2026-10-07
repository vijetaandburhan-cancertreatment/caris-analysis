"""Technical completion inventory, not clinical interpretation or fusion validation."""
import pathlib,json,csv,collections,hashlib,datetime,pysam,re
P=pathlib.Path(__file__).resolve().parent;D=P/'patient-D8'
run=json.loads((D/'run.json').read_text());assert run['status']=='complete'
assert int(run['STAR_final']['Number of input reads'])==23209264
def read_tsv(name):
 lines=(D/name).read_text().splitlines()
 return list(csv.DictReader([lines[0].lstrip('#')]+lines[1:],delimiter='\t'))
accepted=read_tsv('fusions.tsv');discarded=read_tsv('fusions.discarded.tsv')
compact=[]
for n,r in enumerate(accepted,1):
 q={k:r[k] for k in ['gene1','gene2','breakpoint1','breakpoint2','confidence','reading_frame','type','filters','tags','split_reads1','split_reads2','discordant_mates']}
 ids=[] if r['read_identifiers'] in ['.',''] else r['read_identifiers'].split(',')
 q.update(candidate_number=n,reported_unique_read_identifiers=len(set(ids)),reported_split_plus_discordant=sum(int(r[k]) for k in ['split_reads1','split_reads2','discordant_mates']))
 compact.append(q)
target=D/'target-genes.bam';pysam.samtools.quickcheck('-v',str(target))
target_counts=collections.Counter();target_flags=collections.Counter()
with pysam.AlignmentFile(target,'rb') as bam:
 for r in bam:
  target_counts[r.reference_name]+=1
  for key,condition in [('primary',not(r.is_secondary or r.is_supplementary)),('secondary',r.is_secondary),('supplementary',r.is_supplementary),('duplicate',r.is_duplicate)]:
   if condition:target_flags[key]+=1
sorted_target=D/'target-genes.sorted.bam'
if not sorted_target.exists():pysam.samtools.sort('-@','2','-m','256M','-o',str(sorted_target),str(target))
if not pathlib.Path(str(sorted_target)+'.bai').exists():pysam.samtools.index(str(sorted_target))
pysam.samtools.quickcheck('-v',str(sorted_target))
meta={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'complete','input_pairs':23209264,'STAR_final':run['STAR_final'],'accepted_rows':len(accepted),'discarded_rows':len(discarded),'accepted_confidence_counts':dict(collections.Counter(r['confidence'] for r in accepted)),'accepted_inventory':compact,'target_BAM':{'records_by_contig':dict(target_counts),'flags':dict(target_flags),'quickcheck_pass':True,'sorted_indexed_copy':str(sorted_target),'selection':'Overlaps four full gene intervals, including all record flags. Mates outside intervals may be absent; downstream filters still needed.','duplicate_flag_caveat':'Fresh STAR output has not undergone BAM duplicate marking; flag absence does not establish independent molecules.'},'limits':['Caller-reported support and confidence are not direct independent molecular validation.','The sparse D8 index has observed control alignment differences from D1.','Arriba does not report intragenic deletions.','No full patient BAM stored locally.','No matched normal; no clinical validation.'],'files':[]}
log_text=(D/'Arriba.log').read_text();warning=re.search(r'(\d+) SAM records were malformed and ignored',log_text)
meta['caller_warnings']={'malformed_warning_counter':int(warning.group(1)) if warning else 0,'malformed_counter_unit':'Mixed record-level wrong-end clipping events and rejected query-name+HI alignment groups; not a unique-read fraction. See separate diagnostic audit.','subsampling_warning_present':'some fusions were subsampled' in log_text,'subsampling_threshold_default':300,'interpretation_status':'Technical completion only. Warning categories and candidate evidence must be audited before biological conclusions.'}
for f in sorted(D.iterdir()):
 if f.is_file() and f.name not in ['completion-summary.json','candidate-inventory.tsv']:
  with f.open('rb') as h: digest=hashlib.file_digest(h,'sha256').hexdigest()
  meta['files'].append({'name':f.name,'bytes':f.stat().st_size,'sha256':digest})
(D/'completion-summary.json').write_text(json.dumps(meta,indent=2))
if compact:
 with (D/'candidate-inventory.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(compact[0]),delimiter='\t');w.writeheader();w.writerows(compact)
print(json.dumps({k:meta[k] for k in ['status','input_pairs','accepted_rows','discarded_rows','accepted_confidence_counts','target_BAM']},indent=2))

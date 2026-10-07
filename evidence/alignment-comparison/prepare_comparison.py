from pathlib import Path
import json,csv,collections,hashlib,datetime,bisect,pysam,shutil,os
R=Path(__file__).resolve().parent;B=R.parent;P=B/'oct5-population-panel-concordance-v1';W=Path('/Users/burhanazeem/Documents/Codex/2026-09-05/finances-plugin-finances-openai-curated-remote-3')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
cohort=rows(P/'joint-callable.tsv');assert len(cohort)==4015
exceptions=rows(P/'independent-count-audit/all13-below90-dominant-sites.tsv')+rows(P/'independent-count-audit/three-MAPQ60-only-flags.tsv');assert len(exceptions)==16
assert len({(x['chrom'],x['pos1']) for x in exceptions})==16
shutil.copy2(P/'joint-callable.tsv',R/'baseline-fixed-4015.tsv')
with (R/'baseline-sites.bed').open('w') as f:
 for x in cohort:f.write(f"{x['chrom']}\t{int(x['pos1'])-1}\t{x['pos1']}\n")
with (R/'frozen-exception-sites.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(exceptions[0]),delimiter='\t');w.writeheader();w.writerows(exceptions)
names=set();by={};flag=collections.Counter();namecount=collections.Counter()
with pysam.AlignmentFile(str(B/'TN26-279853/RNA_TN26-279853.bam'),'rb',index_filename=str(W/'work/oct1-analysis/RNA_TN26-279853.bam.bai')) as bam:
 header=bam.header.to_dict()
 for x in exceptions:
  pos=int(x['pos1'])-1;site=(x['chrom'],pos);site_names=set();records=0
  for rd in bam.fetch(site[0],pos,pos+1):
   if rd.is_unmapped:continue
   if not any(s<=pos<e for s,e in rd.get_blocks()):continue
   names.add(rd.query_name);site_names.add(rd.query_name);records+=1;flag[str(rd.flag)]+=1;namecount[rd.query_name]+=1
  by[f'{site[0]}:{pos+1}']={'QNAMEs':len(site_names),'Mblock_overlapping_records':records,'capture_selection':'Any aligned M/=/X base covers frozen exception position; no flag/MAPQ/BQ filters. All QNAME records captured in new stream regardless of locus/flag.'}
(R/'frozen-exception-qnames.txt').write_text('\n'.join(sorted(names))+'\n')
(R/'original-header.json').write_text(json.dumps(header,indent=2)+'\n')
summary={'qnames':len(names),'sites':16,'by_site':by,'capture_name_selection_flag_counts':dict(flag),'read_id_handling':'Exact BAM QNAME; no normalization or suffix alteration. Original/new pipeline must show capture name coverage.','input_sha256':{'cohort':sha(R/'baseline-fixed-4015.tsv'),'exception_sites':sha(R/'frozen-exception-sites.tsv'),'exception_qnames':sha(R/'frozen-exception-qnames.txt')},'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'frozen-selection.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

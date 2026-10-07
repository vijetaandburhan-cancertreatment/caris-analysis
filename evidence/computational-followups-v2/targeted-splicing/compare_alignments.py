#!/usr/bin/env python3
import csv,json,pathlib,datetime
BASE=pathlib.Path(__file__).resolve().parent
A=BASE/'original-Caris';B=BASE/'new-D8';O=BASE/'comparison';O.mkdir(exist_ok=True)
def load(p):
 s=json.load(open(p/'summary.json'))
 j={ (r['gene'],int(r['intron_start0']),int(r['intron_end0'])):r for r in csv.DictReader(open(p/'junctions.tsv'),delimiter='\t')}
 f={ (r['gene'],r['fragment'].split('|',1)[1]):r['state'] for r in json.load(open(p/'allele-fragment-evidence.json'))}
 return s,j,f
a,aj,af=load(A);b,bj,bf=load(B)
assert a['filters']==b['filters']
rows=[]
for k in sorted(set(aj)|set(bj)):
 ar=aj.get(k,{});br=bj.get(k,{})
 r=dict(gene=k[0],intron_start0=k[1],intron_end0=k[2],annotation=ar.get('reference_annotation_transcripts',br.get('reference_annotation_transcripts','')))
 for n,src in [('old',ar),('new',br)]:
  for v in ['query_name_fragments','distinct_fragment_alignment_patterns','alternate_linked_pairs','alternate_and_junction_same_read']:r[n+'_'+v]=int(src.get(v,0))
 rows.append(r)
with (O/'junction-comparison.tsv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
genes={}
for g in ('BAP1','RASA1','LATS1','LATS2'):
 aa=next(s for s in a['summaries'] if s['gene']==g);bb=next(s for s in b['summaries'] if s['gene']==g)
 old={k[1]:v for k,v in af.items() if k[0]==g};new={k[1]:v for k,v in bf.items() if k[0]==g}
 common=set(old)&set(new)
 genes[g]={'old':aa,'new':bb,'informative_names_shared':len(common),'informative_names_same_state':sum(old[k]==new[k] for k in common),'old_only':len(set(old)-set(new)),'new_only':len(set(new)-set(old)),'different_state':[{'query_name':k,'old':old[k],'new':new[k]} for k in sorted(common) if old[k]!=new[k]]}
focus={'BAP1 lower mutation-exon flank':('BAP1',52408077,52408473),'BAP1 upper mutation-exon flank':('BAP1',52408606,52409553),'BAP1 mutation-exon skipping':('BAP1',52408077,52409553),'LATS1 lower mutation-exon flank':('LATS1',149684592,149695073),'LATS1 upper mutation-exon flank':('LATS1',149695221,149701778)}
keyed={(r['gene'],r['intron_start0'],r['intron_end0']):r for r in rows}
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'filters':a['filters'],'genes':genes,'focus_junctions':{label:keyed.get(key,{'absent_both':True}) for label,key in focus.items()},'known_limits':['Same specimen and underlying read pool; re-alignment agreement is not independent biological replication.','Raw read-name matching ignores differences in original/new read group labels for this single library.','Junction counts are descriptive query-name counts, not calibrated isoform fractions or independent molecules.','LATS2 has no variant assay in this bounded script; no conclusion about its alternate sequence follows from uncovered state.','BAP1 exon skipping junction is already annotated in GENCODE37 ENST00000490917.1 (nonsense_mediated_decay); it is not novel or proven tumor-specific.']}
(O/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

#!/usr/bin/env python3
"""Inventory all calls and predeclare bounded read-review triage; does not establish actionability."""
import argparse,collections,csv,hashlib,json,pathlib,re
# Genes used only to prioritize review, not to call a particular fusion actionable.
PRIORITY={'ALK','ROS1','RET','NTRK1','NTRK2','NTRK3','NRG1','FGFR1','FGFR2','FGFR3','BRAF','RAF1','MET','NUTM1','SS18','EWSR1','FUS','WWTR1','CAMTA1','TFE3','YAP1','BAP1','LATS1','LATS2','RASA1'}
def read(p,kind):
 with open(p) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   r['gene1']=r.pop('#gene1',r.get('gene1',''))
   r['source_set']=kind
   for field in ['split_reads1','split_reads2','discordant_mates','coverage1','coverage2']:
    r[field]=int(r[field])
   r['split_support']=r['split_reads1']+r['split_reads2'];r['total_reported_support']=r['split_support']+r['discordant_mates']
   r['read_names']=sorted(set(x for x in r.get('read_identifiers','').split(',') if x and x!='.'))
   r['distinct_named_reads']=len(r['read_names'])
   keys=[r.get(k,'') for k in ['gene1','gene2','breakpoint1','breakpoint2','direction1','direction2']]
   r['candidate_id']=hashlib.sha256('|'.join(keys).encode()).hexdigest()[:12]
   genes=set(re.split('[,()]',r['gene1']+','+r['gene2']));r['priority_gene']=sorted(genes&PRIORITY)
   r['review_reasons']=[]
   if kind=='accepted':r['review_reasons'].append('all accepted calls')
   if r['priority_gene'] and r['split_support']>=2:r['review_reasons'].append('priority gene and >=2 split reads')
   if r.get('tags','.')!='.' and r['total_reported_support']>=3:r['review_reasons'].append('caller-supplied database tag and >=3 reported supports')
   if r.get('reading_frame')=='in-frame' and r['split_support']>=3 and r['distinct_named_reads']>=3:r['review_reasons'].append('in-frame with >=3 split and named reads')
   ft=r.get('fusion_transcript','');r['junction_patterns']=[]
   if ft.count('|')==1:
    left,right=ft.split('|');lm=re.search('[ACGTacgt]+$',left);rm=re.match('[ACGTacgt]+',right)
    if lm and rm:
     for arm in [15,20,25]:
      if min(len(lm[0]),len(rm[0]))>=arm:r['junction_patterns'].append({'arm_bases':arm,'sequence':(lm[0][-arm:]+rm[0][:arm]).upper()})
   yield r

def main():
 p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--out',required=True);a=p.parse_args();run=pathlib.Path(a.run);out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
 rows=[]
 for filename,kind in [('fusions.tsv','accepted'),('fusions.discarded.tsv','discarded')]:rows.extend(read(run/filename,kind))
 selected=[r for r in rows if r['review_reasons']]
 (out/'all-candidates.json').write_text(json.dumps(rows,indent=2)+'\n');(out/'review-candidates.json').write_text(json.dumps(selected,indent=2)+'\n')
 fields=['candidate_id','source_set','gene1','gene2','breakpoint1','breakpoint2','type','confidence','reading_frame','split_reads1','split_reads2','discordant_mates','distinct_named_reads','filters','tags','priority_gene','review_reasons']
 with (out/'candidate-triage.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
  for r in rows:w.writerow({**r,'priority_gene':','.join(r['priority_gene']),'review_reasons':'; '.join(r['review_reasons'])})
 stats={'total_rows':len(rows),'accepted':sum(r['source_set']=='accepted' for r in rows),'discarded':sum(r['source_set']=='discarded' for r in rows),'selected_for_read_review':len(selected),'accepted_confidence':dict(collections.Counter(r['confidence'] for r in rows if r['source_set']=='accepted')),'predeclared_triage':'All accepted; discarded priority-gene with >=2split OR database-tagged >=3supports OR in-frame >=3split andnamed. Database tags are not proof. Full discarded output preserved.','limitations':['No actionability classification from gene name alone.','Caller read names/support counts are not independent molecule counts.','Read-level and reference-ambiguity audit must follow; this script only inventories.','Discarded candidates outside bounded triage are not individually validated; a negative review is not comprehensive exclusion.']}
 (out/'inventory.json').write_text(json.dumps(stats,indent=2)+'\n');print(json.dumps(stats,indent=2))
if __name__=='__main__':main()

"""Exact central40-nt join support, allowing differences outside that window."""
import json,pathlib,sys
from Bio import SeqIO
p=pathlib.Path(sys.argv[1]).resolve();pref=str(p)[:-5];d=json.loads(p.read_text());f={s.id:str(s.seq).upper()for s in SeqIO.parse(pref+'.fusions.fasta','fasta')};tab=str.maketrans('ACGTN','TGCAN');rows=[];details=[]
for g in d['genes']:
 names=set();seqs=set();all_ev=[]
 for t in g.get('transcripts',[]):
  a=t['transcriptA'];bp=a['endPos']-a['startPos'];fs=f.get(t['fasta_record'],'')
  if bp<20 or len(fs)-bp<20:continue
  win=fs[bp-20:bp+20]
  if not set(win)<=set('ACGT'):continue
  ev=[]
  for z in g.get('readpairs',[]):
   for mate in ['read1','read2']:
    s=z[mate]['seq'].upper()
    for ori,query in [('forward',s),('reverse-complement',s.translate(tab)[::-1])]:
     start=query.find(win)
     if start<0:continue
     names.add(z['read1']['name']);seqs.add((z['read1']['seq'],z['read2']['seq']));ev.append({'name':z['read1']['name'],'mate':mate,'orientation':ori,'join_offset_in_oriented_read':start+20})
  if ev:all_ev.append({'fasta_record':t['fasta_record'],'read_evidence':ev})
 row={'geneA':g['geneA']['name'],'geneB':g['geneB']['name'],'pizzly_pairs':g['paircount'],'pizzly_splits':g['splitcount'],'junction40_name_pairs':len(names),'junction40_sequence_pairs':len(seqs),'supported_transcript_combinations':len(all_ev)};rows.append(row);details.append(dict(row,junctions=all_ev))
rows.sort(key=lambda x:(-x['junction40_sequence_pairs'],x['geneA'],x['geneB']))
res={'method':'Search a 40nt sequence centered on each proposed join exactly within either read/orientation, preserving20nt each side. Unlike whole-read exact validation, differences elsewhere in the read do not exclude support.','caveats':['Still transcriptome-only and dependent on proposed breakpoint; absence is not a fusion-negative conclusion.','Distinct names/sequences are not independent molecules. Paralogy, normal read-through and FFPE/library chimeras require separate review.'],'summary':rows,'details':details};pathlib.Path(pref+'.junction-window-validation.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({'candidates':len(rows),'top25':rows[:25]},indent=2))

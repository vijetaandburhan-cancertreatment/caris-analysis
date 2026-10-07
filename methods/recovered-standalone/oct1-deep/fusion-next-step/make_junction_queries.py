import json,pathlib,sys
from Bio import SeqIO
p=pathlib.Path(sys.argv[1]);d=json.loads(p.read_text());f={r.id:str(r.seq).upper() for r in SeqIO.parse(str(p)[:-5]+'.fusions.fasta','fasta')}
tr=str.maketrans('ACGT','TGCA');queries=[];metadata=[]
for gi,g in enumerate(d['genes']):
 for ti,t in enumerate(g.get('transcripts',[])):
  seq=f.get(t['fasta_record'],'');a=t['transcriptA'];bp=a['endPos']-a['startPos']
  if bp<20 or len(seq)-bp<20:continue
  s=seq[bp-20:bp+20]
  if not set(s)<=set('ACGT'):continue
  qid=f'gene{gi}_transcript{ti}'
  for orient,q in [('fwd',s),('rev',s.translate(tr)[::-1])]:queries.append((qid+'_'+orient,q))
  metadata.append({'query_id':qid,'geneA':g['geneA']['name'],'geneB':g['geneB']['name'],'fasta_record':t['fasta_record'],'join0':bp})
prefix=str(p)[:-5]+'.junction40'
pathlib.Path(prefix+'.queries.tsv').write_text(''.join(a+'\t'+b+'\n' for a,b in queries))
pathlib.Path(prefix+'.metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps({'junctions':len(metadata),'orientation_queries':len(queries),'prefix':prefix}))

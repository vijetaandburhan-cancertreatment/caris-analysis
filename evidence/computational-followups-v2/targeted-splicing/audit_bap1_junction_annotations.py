import pathlib,re,csv,json,collections,pysam
W=pathlib.Path(__file__).resolve().parent;C=pathlib.Path.home()/'.local/share/codex/caris-analysis/genome-fusion-reference';gtf=C/'gencode.v37.primary_assembly.annotation.gtf'
cds=[];exons=collections.defaultdict(list);types={};mane=[]
with gtf.open() as f:
 for line in f:
  if 'gene_name "BAP1"' not in line:continue
  t=line.rstrip().split('\t');a=dict(re.findall(r'(\w+) "([^"]*)"',t[8]));tx=a.get('transcript_id');start,end=int(t[3])-1,int(t[4])
  if tx and a.get('transcript_type'):types[tx]=a['transcript_type']
  if t[2]=='transcript' and 'tag "MANE_Select"' in t[8]:mane.append(tx)
  if t[2]=='CDS' and tx=='ENST00000460680.6':cds.append((start,end))
  if t[2]=='exon':exons[tx].append((start,end))
known=collections.defaultdict(list)
for tx,intervals in exons.items():
 z=sorted(set(intervals))
 for a,b in zip(z,z[1:]):known[a[1],b[0]].append(tx)
fa=pysam.FastaFile(str(C/'GRCh38.primary_assembly.genome.fa'))
rc=lambda x:x.translate(str.maketrans('ACGT','TGCA'))[::-1]
rows={}
for label in ['original-Caris','new-D8']:
 rows[label]={(int(r['intron_start0']),int(r['intron_end0'])):r for r in csv.DictReader((W/label/'junctions.tsv').open(),delimiter='\t') if r['gene']=='BAP1'}
report={'MANE_Select_transcripts':mane,'MANE_reference_used':'ENST00000460680.6','CDS_intervals0':sorted(cds),'junctions':[]}
for a,b in [(52402871,52403182),(52404586,52405764),(52408077,52409553),(52408077,52408473),(52408606,52409553)]:
 left=fa.fetch('chr3',a,a+2).upper();right=fa.fetch('chr3',b-2,b).upper();n=sum(max(0,min(b,y)-max(a,x)) for x,y in cds)
 entry={'intron_start0':a,'intron_end0':b,'forward_motif':left+'-'+right,'BAP1_strand_donor_acceptor':rc(right)+'-'+rc(left),'MANE_CDS_overlap_bases':n,'CDS_overlap_modulo3':n%3,'annotation':sorted(known[a,b]),'annotation_types':{tx:types[tx] for tx in known[a,b]}}
 for label in rows:
  r=rows[label][a,b];entry[label]={k:int(r[k]) for k in ['query_name_fragments','alternate_linked_pairs','alternate_and_junction_same_read']}
 report['junctions'].append(entry)
assert mane==['ENST00000460680.6']
assert [e['MANE_CDS_overlap_bases'] for e in report['junctions'][:3]]==[45,185,133]
assert all(e['BAP1_strand_donor_acceptor']=='GT-AG' for e in report['junctions'])
assert not report['junctions'][0]['annotation'] and not report['junctions'][1]['annotation']
assert report['junctions'][2]['annotation_types']=={'ENST00000490917.1':'nonsense_mediated_decay'}
(W/'independent-BAP1-junction-annotations.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

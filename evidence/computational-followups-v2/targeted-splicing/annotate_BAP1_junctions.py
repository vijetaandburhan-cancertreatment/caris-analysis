import csv,json,pathlib,re
import pysam
B=pathlib.Path(__file__).resolve().parent
R=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference')
fa=pysam.FastaFile(str(R/'GRCh38.primary_assembly.genome.fa'))
cds=[];txmeta={}
for line in open(R/'gencode.v37.primary_assembly.annotation.gtf'):
 if 'gene_name "BAP1"' not in line:continue
 a=line.strip().split('\t');d=dict(re.findall(r'(\w+) "([^"]*)"',a[8]))
 if a[2]=='transcript':txmeta[d['transcript_id']]={k:d.get(k) for k in ('transcript_name','transcript_type','transcript_support_level')}
 if a[2]=='CDS' and d.get('transcript_id')=='ENST00000460680.6':cds.append((int(a[3])-1,int(a[4])))
rows=[]
for r in csv.DictReader(open(B/'comparison/junction-comparison.tsv'),delimiter='\t'):
 if r['gene']!='BAP1':continue
 s,e=int(r['intron_start0']),int(r['intron_end0']);left=fa.fetch('chr3',s,s+2).upper();right=fa.fetch('chr3',e-2,e).upper()
 # BAP1 is minus strand; donor/acceptor in transcription direction.
 comp=str.maketrans('ACGT','TGCA');rc=lambda x:x.translate(comp)[::-1]
 lost=sum(max(0,min(e,b)-max(s,a)) for a,b in cds)
 r.update(forward_reference_splice_edges=left+'-'+right,BAP1_strand_donor_acceptor=rc(right)+'-'+rc(left),MANE_CDS_bases_inside_this_intron=lost,MANE_CDS_removed_modulo3=lost%3,annotation_types=';'.join(t+':'+str(txmeta.get(t,{}).get('transcript_type')) for t in r['annotation'].split(',') if t))
 rows.append(r)
with open(B/'comparison/BAP1-junction-annotations.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps([r for r in rows if not r['annotation'] and int(r['old_query_name_fragments'])>=30 and int(r['new_query_name_fragments'])>=30],indent=2))
(B/'comparison/BAP1-junction-annotation-limits.txt').write_text('GENCODE37 MANE Select ENST00000460680.6. Coding bases inside a proposed intron are a sequence-level description relative to this reference, not a reconstruction of a full-length RNA/protein. Events can be physiological, cell-mixture-derived, technically artifactual or belong to another isoform; lack of G37 annotation does not make an event cancer-specific. The downstream junctions are not phased to the early BAP1 frameshift, and cannot establish functional rescue, NMD escape, an antigen or drug sensitivity. Counts are query-name fragments, not independent molecules.\n')

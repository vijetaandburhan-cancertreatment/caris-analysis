"""Independent exact-string validation of Pizzly local sequences; no fusion calling.

This intentionally does not use Pizzly's transcript read-index arrays. Each actual
read sequence is independently searched (both orientations) against each reported
fused transcript sequence. Exact whole-read matches must span the reported join
with >=20 nt on each side. These conservative positive checks have no sensitivity
claim; a failed exact match is not evidence of absence.
"""
import csv, gzip, json, pathlib, sys
from Bio import SeqIO

source = pathlib.Path(sys.argv[1]).resolve()
root = pathlib.Path(__file__).resolve().parent
prefix = str(source)[:-5]
fused = {r.id:str(r.seq).upper() for r in SeqIO.parse(prefix+'.fusions.fasta','fasta')}
d = json.loads(source.read_text())
needed = {t[p]['id'] for g in d['genes'] for t in g.get('transcripts',[]) for p in ['transcriptA','transcriptB']}
normal = {}
with gzip.open(root/'public/gencode.v37.transcripts.fa.gz','rt') as f:
    for r in SeqIO.parse(f,'fasta'):
        rid=r.id.split('|')[0]
        if rid in needed: normal[rid]=str(r.seq).upper()
tr = str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(tr)[::-1]
def allfind(seq,query):
    at=seq.find(query)
    while at>=0:
        yield at
        at=seq.find(query,at+1)
rows=[]; details=[]
for g in d['genes']:
    confirmed_pairs=set();confirmed_sequences=set();normal_matches=set();observed=set();txdetails=[]
    for t in g.get('transcripts',[]):
        fname=t['fasta_record']; fs=fused.get(fname)
        if fs is None: continue
        ta,tb=t['transcriptA'],t['transcriptB']
        bp=ta['endPos']-ta['startPos']
        if bp<=0 or bp>=len(fs):continue
        txnames=set();txseqs=set(); starts=set(); normalnames=set(); evidence=[]
        partners=[normal.get(ta['id'],''),normal.get(tb['id'],'')]
        for rp in g.get('readpairs',[]):
            key=(rp['read1']['name'],rp['read2']['name']); sq=(rp['read1']['seq'],rp['read2']['seq'])
            for mate in ['read1','read2']:
                read=rp[mate];seq=read['seq'].upper()
                for orient,s in [('forward',seq),('reverse-complement',rc(seq))]:
                    for start in allfind(fs,s):
                        left=bp-start;right=start+len(s)-bp
                        if left<20 or right<20:continue
                        txnames.add(key);txseqs.add(sq);starts.add((start,orient,mate))
                        same_normal=any(s in p or rc(s) in p for p in partners if p)
                        if same_normal:normalnames.add(key)
                        evidence.append({'pair_name':key,'mate':mate,'orientation':orient,'start0':start,'left_anchor':left,'right_anchor':right,'whole_read_matches_one_normal_parent':same_normal})
        if txnames:
            confirmed_pairs |= txnames; confirmed_sequences |= txseqs; normal_matches |=normalnames
            observed.add(fname)
            txdetails.append({'fasta_record':fname,'join0':bp,'exact_spanning_name_pairs':len(txnames),'exact_spanning_sequence_pairs':len(txseqs),'distinct_start_orientation_mate':len(starts),'normal_parent_match_name_pairs':len(normalnames),'evidence':evidence})
    row={'geneA':g['geneA']['name'],'geneB':g['geneB']['name'],'pizzly_pairs':g['paircount'],'pizzly_splits':g['splitcount'],'exact20nt_join_name_pairs':len(confirmed_pairs),'exact20nt_join_sequence_pairs':len(confirmed_sequences),'name_pairs_also_matching_normal_parent':len(normal_matches),'exact_supported_isoform_combinations':len(observed)}
    rows.append(row);details.append(dict(row,transcripts=txdetails))
rows.sort(key=lambda x:(-x['exact20nt_join_sequence_pairs'],x['geneA'],x['geneB']))
notes=['Whole-read exact match across Pizzly-reported join with >=20 nt each side; independent string search of reported sequences, not transcript read-index arrays.','Negative exact-check result is not a fusion-negative result; mismatch, FFPE damage or imperfect breakpoint may prevent an exact match.','Normal-parent check covers the two reported partner transcripts, not all isoforms, paralogs or the genome.','No UMI tags or RNA duplicate flags: distinct names and sequence pairs are not independent molecule counts.','Unvalidated research candidates require orthogonal junction testing, genomic mapping and technical/normal-readthrough review.']
out=pathlib.Path(prefix+'.exact-junctions')
pathlib.Path(str(out)+'.json').write_text(json.dumps({'source':str(source),'minimum_each_side_anchor_nt':20,'notes':notes,'summary':rows,'candidate_details':details},indent=2)+'\n')
with open(str(out)+'.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['geneA','geneB'],delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({'candidate_count':len(rows),'top25':rows[:25],'notes':notes},indent=2))

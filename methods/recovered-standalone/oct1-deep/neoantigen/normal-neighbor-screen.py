"""BAP1 peptide sequence neighbors in the human reference proteome.

Not an off-target or T-cell cross-reactivity assay. Complete search only for
zero, one or two substitutions at equal peptide length in this reference.
"""
from pathlib import Path
from collections import defaultdict
import gzip,json,time,csv
from Bio import SeqIO
OUT=Path(__file__).resolve().parent
QUERIES=['IEERKGLYL','EERKGLYL']
started=time.time();hits=defaultdict(list);proteins=0
with gzip.open(OUT/'gencode.v37.pc_translations.fa.gz','rt') as f:
    for record in SeqIO.parse(f,'fasta'):
        proteins+=1;seq=str(record.seq)
        for query in QUERIES:
            n=len(query);cuts=[0,n//3,2*n//3,n]
            starts=set()
            # At most two substitutions leave at least one of these three
            # disjoint seed segments unchanged (pigeonhole principle).
            for i in range(3):
                offset=cuts[i];seed=query[offset:cuts[i+1]];at=seq.find(seed)
                while at>=0:
                    s=at-offset
                    if 0<=s<=len(seq)-n:starts.add(s)
                    at=seq.find(seed,at+1)
            for s in starts:
                peptide=seq[s:s+n];diff=[i+1 for i,(a,b) in enumerate(zip(query,peptide)) if a!=b]
                if len(diff)<=2:
                    hits[(query,peptide)].append({'reference_header':record.description,'start_aa1':s+1,'end_aa1':s+n,'mismatch_positions1':diff,'left_flank5':seq[max(0,s-5):s],'right_flank5':seq[s+n:s+n+5]})
rows=[]
for (query,peptide),where in sorted(hits.items()):
    rows.append({'query':query,'normal_reference_peptide':peptide,'substitutions':len(where[0]['mismatch_positions1']),'mismatch_positions1':where[0]['mismatch_positions1'],'reference_records':len(where),'occurrences':where})
(OUT/'bap1-normal-sequence-neighbors.json').write_text(json.dumps({'reference':'GENCODE37 protein FASTA, same independently hashed file as main analysis','reference_records_scanned':proteins,'scope':'Equal-length exact/one-/two-substitution neighbors of the two BAP1 candidates; all candidate seeds searched; no gaps, no normal population polymorphisms, no post-translational/noncanonical ORFs, no TCR recognition or safety inference.','elapsed_seconds':time.time()-started,'results':rows},indent=2)+'\n')
with (OUT/'bap1-normal-sequence-neighbors.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['BAP1_candidate','normal_reference_neighbor','substitutions','mismatch_positions1','reference_record_count','example_reference_header'])
    for r in rows:w.writerow([r['query'],r['normal_reference_peptide'],r['substitutions'],','.join(map(str,r['mismatch_positions1'])),r['reference_records'],r['occurrences'][0]['reference_header']])
print('Scanned',proteins,'proteins; unique query/neighbor pairs',len(rows),'seconds',time.time()-started)
for r in rows:print(r['query'],r['normal_reference_peptide'],r['substitutions'],r['reference_records'],r['occurrences'][0]['reference_header'])

"""Describe premature-stop/exon geometry; no NMD probability or escape call."""
from pathlib import Path
import json,gzip,re
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
base=ROOT/'work/oct1-deep/coding-indels'
rows=[r for r in json.loads((base/'conditional-transcript-consequences.json').read_text())['results'] if r['gene'] in ['BAP1','RASA1','LATS1']]
ann={r['transcript']:r for r in json.loads((base/'selected-coding-transcripts.json').read_text())}
exons={r['transcript']:[] for r in rows}
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'):continue
        a=line.rstrip().split('\t')
        if a[2]!='exon':continue
        m=re.search(r'transcript_id "([^"]+)"',a[8])
        if m and m[1] in exons:exons[m[1]].append((int(a[3])-1,int(a[4])))
results=[]
for r in rows:
    tx=r['transcript'];a=ann[tx];ex=sorted(exons[tx],reverse=a['strand']=='-')
    mrna=[p for s,e in ex for p in (range(s,e) if a['strand']=='+' else range(e-1,s-1,-1))]
    coding=[p for s,e in sorted(a['CDS'],reverse=a['strand']=='-') for p in (range(s,e) if a['strand']=='+' else range(e-1,s-1,-1))]
    offset=mrna.index(coding[0]);edit_start=offset+r['coding_edit_start0'];edit_end=offset+r['coding_edit_end0'];delta=len(r['coding_alt'])-len(r['coding_ref'])
    # Map each mutant nucleotide separately: the LATS1 stop includes an
    # inserted base, which has no reference coordinate.
    stop_start=offset+r['mutant_protein_length']*3;stop_end=stop_start+3
    mutant_map=mrna[:edit_start]+[None]*len(r['coding_alt'])+mrna[edit_end:]
    stop_positions=mutant_map[stop_start:stop_end]
    assert len(stop_positions)==3
    bounds=[];cumulative=0
    for s,e in ex[:-1]:
        cumulative+=e-s
        assert not edit_start<cumulative<edit_end
        bounds.append(cumulative+(delta if cumulative>=edit_end else 0))
    downstream=[b for b in bounds if b>=stop_end]
    genomic_stop_base=next(p for p in stop_positions if p is not None)
    exon_index=next(i+1 for i,(s,e) in enumerate(ex) if s<=genomic_stop_base<e)
    result={'gene':r['gene'],'transcript':tx,'exon_count_including_UTRs':len(ex),'stop_exon_number_transcript_order':exon_index,'conditional_stop_codon_mutant_CDS_1based':[r['mutant_protein_length']*3+1,r['mutant_protein_length']*3+3],'conditional_stop_reference_CDS_1based':[mrna.index(p)-offset+1 if p is not None else None for p in stop_positions],'stop_codon_genome_positions1':[p+1 if p is not None else None for p in stop_positions],'distance_stop_end_to_next_downstream_junction_nt':min(downstream)-stop_end if downstream else None,'distance_stop_end_to_final_junction_nt':bounds[-1]-stop_end,'downstream_junction_count':len(downstream),'new_tail_amino_acids':r['new_tail_length'],'interpretation':'Descriptive reference-isoform geometry only. Null reference coordinate denotes inserted base. No categorical NMD-escape claim, no NMD efficiency probability, no full-length protein/presentation claim. Actual mutant RNA is observed under separate read analyses.'}
    results.append(result);print(json.dumps(result),flush=True)
(OUT/'stop-context.json').write_text(json.dumps({'scope':'Conditional single-edit reference initiation, GENCODE37 selected transcript, splice boundaries shifted by exact indel length. Stop start/end checked against edit. No use of variant position as surrogate for actual new stop.','results':results},indent=2)+'\n')

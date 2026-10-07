"""Conditional single-allele GENCODE37 transcript translation.

This intentionally does not infer a complete expressed haplotype. Genomic
indel length alone is insufficient near exon boundaries or compound variants.
"""
from pathlib import Path
from collections import defaultdict
import gzip,json,re,csv
from Bio.Seq import Seq
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
CAND=json.loads((OUT/'expressed-frameshift-shortlist.json').read_text())['candidates']
ANN={r['transcript']:r for r in json.loads((OUT/'selected-coding-transcripts.json').read_text())}
RNA=json.loads((OUT/'selected-reference-transcripts.json').read_text())
PROTEIN=json.loads((OUT/'selected-reference-proteins.json').read_text())
EXONS=defaultdict(list)
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'):continue
        a=line.rstrip('\n').split('\t')
        if len(a)!=9 or a[2]!='exon':continue
        m=re.search(r'transcript_id "([^"]+)"',a[8])
        if m and m[1] in RNA:EXONS[m[1]].append([int(a[3])-1,int(a[4])])

def consequence(v,tx):
    a=ANN[tx];exons=sorted(EXONS[tx],reverse=a['strand']=='-');rnaseq=RNA[tx]['sequence']
    mrnapos=[p for s,e in exons for p in (range(s,e) if a['strand']=='+' else range(e-1,s-1,-1))]
    assert len(mrnapos)==len(rnaseq),(tx,len(mrnapos),len(rnaseq))
    mmap={p:i for i,p in enumerate(mrnapos)}
    cds=sorted(a['CDS'],reverse=a['strand']=='-')
    cpos=[p for s,e in cds for p in (range(s,e) if a['strand']=='+' else range(e-1,s-1,-1))]
    cmap={p:i for i,p in enumerate(cpos)}
    offset=mmap[cpos[0]];coding=rnaseq[offset:offset+len(cpos)]
    assert [mmap[p] for p in cpos]==list(range(offset,offset+len(cpos))),tx
    translated=str(Seq(coding).translate())
    result={'gene':a['gene'],'transcript':tx,'transcript_priority':a['priority'],'strand':a['strand'],'reference_protein_matches_GENCODE':translated==PROTEIN[tx]['protein'],'reference_protein_length':len(translated),'reference_internal_stop':('*' in translated),'pos1':v['pos1'],'chrom':v['chrom'],'ref':v['ref'],'alt':v['alt']}
    p=v['pos1']-1;ref=v['ref'];alt=v['alt']
    while ref and alt and ref[0]==alt[0]:p+=1;ref=ref[1:];alt=alt[1:]
    while ref and alt and ref[-1]==alt[-1]:ref=ref[:-1];alt=alt[:-1]
    if ref:
        poses=list(range(p,p+len(ref)))
        if any(x not in cmap for x in poses):
            result['status']='not_wholly_within_selected_CDS; genomic indel length is not a coding consequence';return result
        indices=sorted(cmap[x] for x in poses);i=indices[0];j=indices[-1]+1
        if indices!=list(range(i,j)):result['status']='noncontiguous_CDS_edit';return result
    else:
        if p-1 not in cmap or p not in cmap or abs(cmap[p-1]-cmap[p])!=1:
            result['status']='insertion_at_CDS_or_splice_boundary; consequence_unresolved';return result
        i=j=min(cmap[p-1],cmap[p])+1
    cref=str(Seq(ref).reverse_complement()) if a['strand']=='-' else ref
    calt=str(Seq(alt).reverse_complement()) if a['strand']=='-' else alt
    if coding[i:j]!=cref:
        result['status']='reference_transcript_genome_mismatch';result['transcript_reference_at_edit']=coding[i:j];result['expected_reference_at_edit']=cref;return result
    mutant_translated_region=rnaseq[offset:offset+i]+calt+rnaseq[offset+j:]
    # Translation is deliberately based on the full 3-prime transcript to stop,
    # retaining UTR-derived codons when the edit changes the terminal frame.
    mutant=str(Seq(mutant_translated_region[:len(mutant_translated_region)//3*3]).translate(to_stop=True))
    prefix=0
    while prefix<min(len(translated),len(mutant)) and translated[prefix]==mutant[prefix]:prefix+=1
    distance=min(min(abs(p-x),abs(p+max(len(ref),1)-1-x)) for s,e in cds for x in (s,e-1))
    result.update(status='conditional_single_edit_translation',coding_edit_start0=i,coding_edit_end0=j,coding_ref=cref,coding_alt=calt,mutant_protein_length=len(mutant),first_changed_or_lost_residue1=prefix+1,reference_context=translated[max(0,prefix-10):prefix+25],mutant_context=mutant[max(0,prefix-10):prefix+40],new_tail_length=len(mutant)-prefix,new_tail_sequence=mutant[prefix:],reference_protein=translated,mutant_protein=mutant,distance_to_nearest_CDS_exon_edge=distance,near_exon_edge=distance<=2)
    return result

if __name__=='__main__':
    rows=[consequence(v,t) for v in CAND for t in v['selected_transcripts']]
    (OUT/'conditional-transcript-consequences.json').write_text(json.dumps({'method':'One GENCODE37 selected transcript per gene; exon/CDS mapping against official GENCODE37 transcript FASTA; reference protein independently checked against official translation FASTA; minimal genomic allele projected into coding orientation. Mutant protein assumes reference initiation and only this single edit, translating through the3-prime transcript to first stop. Not a complete tumor haplotype, isoform or protein/presentation measurement. Boundary edits and compound variation require separate treatment.','results':rows},indent=2)+'\n')
    fields=['gene','transcript','transcript_priority','chrom','pos1','ref','alt','status','reference_protein_matches_GENCODE','reference_protein_length','mutant_protein_length','first_changed_or_lost_residue1','new_tail_length','near_exon_edge','reference_context','mutant_context']
    with (OUT/'conditional-transcript-consequences.tsv').open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader()
        for r in rows:w.writerow({k:r.get(k,'') for k in fields})
    for r in rows:print(r['gene'],r['status'],'aa',r.get('first_changed_or_lost_residue1'),'tail',r.get('new_tail_length'),'edge',r.get('near_exon_edge'),'reference_ok',r['reference_protein_matches_GENCODE'],flush=True)

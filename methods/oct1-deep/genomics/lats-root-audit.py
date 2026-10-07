"""Root independent CIGAR-walking recount and exon-edit protein reconstruction."""
import json,gzip,csv,collections,statistics
from pathlib import Path
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
SRC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
variants=[dict(gene='LATS1',chrom='chr6',pos=149695191,ref='T',alt='TA'),dict(gene='LATS2',chrom='chr13',pos=20975259,ref='C',alt='T')]
annotation=json.load((OUT/'target-annotation.json').open())['genes']
proteins={}
with gzip.open(OUT.parent/'neoantigen/gencode.v37.pc_translations.fa.gz','rt') as f:
    for r in SeqIO.parse(f,'fasta'):
        tx=r.id.split('|')[1]
        if tx in [annotation[g]['selected_transcript'] for g in ['LATS1','LATS2']]:proteins[tx]=str(r.seq)
results=[]
for v in variants:
    ref=json.load((OUT/(v['gene']+'.hg38-reference.json')).open());dna=ref['dna'].upper();rs=ref['start'];p=v['pos']-1
    assert dna[p-rs:p-rs+len(v['ref'])]==v['ref']
    a=annotation[v['gene']];normal=[];edited=[];edits=0
    for st,en in a['selected_CDS']:
        seq=dna[st-rs:en-rs];normal.append(seq)
        if st<=p<en:
            seq=seq[:p-st]+v['alt']+seq[p-st+len(v['ref']):];edits+=1
        edited.append(seq)
    assert edits==1
    wt=''.join(normal);mt=''.join(edited)
    if a['strand']=='-':wt=str(Seq(wt).reverse_complement());mt=str(Seq(mt).reverse_complement())
    wtprot=str(Seq(wt).translate());mtprot=str(Seq(mt[:len(mt)//3*3]).translate(to_stop=True))
    assert wtprot.rstrip('*')==proteins[a['selected_transcript']].rstrip('*')
    first=next((i for i,(x,y) in enumerate(zip(wtprot,mtprot)) if x!=y),min(len(wtprot),len(mtprot)))
    item={**v,'transcript':a['selected_transcript'],'reference_protein_matches_GENCODE37':True,'first_altered_position':first+1,'wt_amino_acid':wtprot[first],'mutant_amino_acid':mtprot[first] if first<len(mtprot) else '*','predicted_mutant_length':len(mtprot),'normal_length':len(wtprot),'normal_context':wtprot[max(0,first-10):first+15],'mutant_context':mtprot[max(0,first-10):first+15],'assays':{}}
    for kind in ['DNA','RNA']:
        st=p-8;en=p+len(v['ref'])+8;normal_hap=dna[st-rs:en-rs];alt_hap=normal_hap[:8]+v['alt']+normal_hap[8+len(v['ref']):]
        calls=collections.Counter();frag=collections.defaultdict(set);strands=collections.Counter();qualities=[]
        with pysam.AlignmentFile(str(SRC/f'{kind}_TN26-279853.bam'),index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai')) as b:
            for r in b.fetch(v['chrom'],st,en):
                if r.flag&(4|256|512|1024|2048) or r.mapping_quality<20 or r.query_qualities is None:continue
                genome=r.reference_start;q=0;letters=[];bq=[];endpoints=set()
                for op,n in r.cigartuples:
                    if op in (0,7,8):
                        left=max(st,genome);right=min(en,genome+n)
                        if right>left:
                            q0=q+left-genome;q1=q+right-genome
                            letters.extend(r.query_sequence[q0:q1]);bq.extend(r.query_qualities[q0:q1])
                            if left==st:endpoints.add(st)
                            if right==en:endpoints.add(en-1)
                        genome+=n;q+=n
                    elif op==1:
                        if st<genome<en:letters.extend(r.query_sequence[q:q+n]);bq.extend(r.query_qualities[q:q+n])
                        q+=n
                    elif op in (2,3):genome+=n
                    elif op==4:q+=n
                    elif op in (5,6):pass
                    else:raise ValueError(op)
                if endpoints!={st,en-1} or not bq or min(bq)<20:continue
                hap=''.join(letters);call='ref' if hap==normal_hap else 'alt' if hap==alt_hap else 'other'
                calls[call]+=1;frag[(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name)].add(call)
                strands[call+('_reverse' if r.is_reverse else '_forward')]+=1
                if call=='alt':qualities.append(min(bq))
        fc=collections.Counter(next(iter(c)) if len(c)==1 else 'discordant' for c in frag.values())
        item['assays'][kind]={'reads':dict(calls),'fragments':dict(fc),'strands':dict(strands),'alt_minimum_full_haplotype_BQ':min(qualities) if qualities else None}
    results.append(item)
original=json.load((OUT/'lats-independent-haplotype-recount.json').open())['results']
for v in results:
    for kind,x in v['assays'].items():
        prior=next(a for a in original if a['gene']==v['gene'] and a['pos1']==v['pos'] and a['kind']==kind)
        assert x['reads']==prior['reads'] and x['fragments']==prior['fragments'],(v['gene'],kind,x,prior)
        x['independent_exact_recount_agrees']=True
report={'method':'Separate CIGAR-walking local-haplotype assembly; MAPQ/BQ>=20; duplicate/secondary/supplementary/QCfail exclusion; grouped by RG+queryname. Independent per-exon genomic edit, strand reversal and translation checked against GENCODE37 protein FASTA. Same data and same public reference are not independent clinical validation.','findings':results,'limitations':['Not somatic/germline classified; no normal specimen.','Possible library/alignment artifacts not excluded by recount.','LATS2 missense functional significance unknown.','No treatment selection or diagnosis established by this screen.']}
(OUT/'lats-root-independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

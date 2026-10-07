"""Summarize Pizzly candidates without promoting calls to validated fusions.

No read sequence or QNAME is put in the compact table. Underlying Pizzly files retain
the local sequence evidence. Distinct sequence pairs are not UMI molecule counts.
"""
import collections, csv, gzip, json, pathlib, re, sys

path = pathlib.Path(sys.argv[1]).resolve()
root = pathlib.Path(__file__).resolve().parent
data = json.loads(path.read_text())
genes = {g[k]['id'] for g in data['genes'] for k in ['geneA','geneB']}
annot = {}
with gzip.open(root.parent/'genomics/gencode.v37.annotation.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        v=line.rstrip('\n').split('\t')
        if len(v)!=9 or v[2]!='gene':
            continue
        a=dict(re.findall(r'(\w+) "([^"]*)"',v[8]))
        if a.get('gene_id') in genes:
            annot[a['gene_id']]={'chrom':v[0],'start1':int(v[3]),'end1':int(v[4]),'strand':v[6],'gene_type':a.get('gene_type','')}
rows=[]
for g in data['genes']:
    a,b=g['geneA'],g['geneB']; ga,gb=annot.get(a['id'],{}),annot.get(b['id'],{})
    rp=g.get('readpairs',[])
    names={(x['read1']['name'],x['read2']['name']) for x in rp}
    seqs={(x['read1']['seq'],x['read2']['seq']) for x in rp}
    split=[x for x in rp if x['type']=='SPLIT']
    split_names={(x['read1']['name'],x['read2']['name']) for x in split}
    split_seqs={(x['read1']['seq'],x['read2']['seq']) for x in split}
    gap=''; topology='different chromosomes or missing annotation'
    if ga and gb and ga['chrom']==gb['chrom']:
        gap=max(ga['start1'],gb['start1'])-min(ga['end1'],gb['end1'])-1
        topology='overlapping genes' if gap<0 else ('same chromosome, same strand' if ga['strand']==gb['strand'] else 'same chromosome, opposite strand')
    caution=[]
    if len(split_seqs)<2: caution.append('fewer than2 distinct split sequence-pairs')
    if len(seqs)<3: caution.append('fewer than3 distinct sequence-pairs')
    if rp and len(seqs)<len(rp)/2:caution.append('most supporting entries repeat identical sequence-pairs')
    if isinstance(gap,int) and gap<200000 and ga.get('strand')==gb.get('strand'):caution.append('nearby same-strand genes: evaluate normal read-through/overlap')
    if 'pseudogene' in ga.get('gene_type','') or 'pseudogene' in gb.get('gene_type',''):caution.append('pseudogene partner: mapping ambiguity')
    rows.append({'geneA':a['name'],'geneB':b['name'],'geneA_id':a['id'],'geneB_id':b['id'],
        'pizzly_paircount':g['paircount'],'pizzly_splitcount':g['splitcount'],
        'readpair_entries':len(rp),'distinct_name_pairs':len(names),'distinct_sequence_pairs':len(seqs),
        'distinct_split_name_pairs':len(split_names),'distinct_split_sequence_pairs':len(split_seqs),
        'candidate_transcript_combinations':len(g.get('transcripts',[])),
        'geneA_chrom':ga.get('chrom',''),'geneB_chrom':gb.get('chrom',''),
        'geneA_strand':ga.get('strand',''),'geneB_strand':gb.get('strand',''),
        'genomic_intergene_gap':gap,'topology':topology,'cautions':'; '.join(caution),
        'status':'Unvalidated transcriptome-only candidate; no clinical or neoantigen claim'})
rows.sort(key=lambda r:(-r['distinct_split_sequence_pairs'],-r['distinct_sequence_pairs'],r['geneA'],r['geneB']))
dest=path.parent/(path.stem+'.summary')
dest.with_suffix(dest.suffix+'.json').write_text(json.dumps({'source':str(path),'candidate_count':len(rows),'notes':['Counts are paired-read and sequence-diversity checks, not independent RNA molecule counts.','Pizzly transcript read-index arrays are not used for breakpoint attribution; actual supporting sequences must be rechecked.','Nearby genes can reflect physiological read-through; none of these candidates are validated.'],'candidates':rows},indent=2)+'\n')
if rows:
    with open(dest.with_suffix(dest.suffix+'.tsv'),'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({'candidate_count':len(rows),'top15':rows[:15]},indent=2))

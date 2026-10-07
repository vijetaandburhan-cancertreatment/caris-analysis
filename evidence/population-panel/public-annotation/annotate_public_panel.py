"""GENCODE37 overlap annotations for a frozen PUBLIC SNP universe; no patient input."""
from pathlib import Path
import bisect, collections, csv, datetime, hashlib, json, re, resource, time

OUT=Path(__file__).resolve().parent
PANEL=OUT.parent/'input/public-common-SNP-panel.tsv'
GTF=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference/gencode.v37.primary_assembly.annotation.gtf')
SOURCE=GTF.with_name(GTF.name+'.gz.manifest.json')

def sha(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

def indices_in_closed_interval(sorted_positions, start1, end1):
    """GTF and VCF are both 1-based, inclusive; no conversion is needed."""
    assert 1 <= start1 <= end1
    return range(bisect.bisect_left(sorted_positions,start1),bisect.bisect_right(sorted_positions,end1))

def boundary_tests():
    pos=[9,10,19,20,21,29,30,31]
    a=[pos[i] for i in indices_in_closed_interval(pos,10,20)]
    b=[pos[i] for i in indices_in_closed_interval(pos,20,30)]
    assert a==[10,19,20] and b==[20,21,29,30]
    got={p:sum(s<=p<=e for s,e in [(10,20),(20,30)]) for p in pos}
    assert got=={9:0,10:1,19:1,20:2,21:1,29:1,30:1,31:0}
    assert [1][next(iter(indices_in_closed_interval([1],1,1)))]==1
    return {'synthetic_closed_interval_boundaries':'PASS','overlapping_gene_endpoint':'PASS','single_base_interval':'PASS','positions':pos,'overlap_counts':got}

def main():
    started=time.time(); checks=boundary_tests()
    with PANEL.open(newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==165782
    loci=[(r['chrom'],int(r['pos1'])) for r in rows]
    assert len(set(loci))==len(loci), 'The frozen input contains duplicate genomic loci'
    bychrom=collections.defaultdict(list)
    for i,(chrom,pos) in enumerate(loci): bychrom[chrom].append((pos,i))
    for chrom in bychrom: bychrom[chrom].sort()
    positions={c:[p for p,i in a] for c,a in bychrom.items()}
    rowids={c:[i for p,i in a] for c,a in bychrom.items()}
    spans=[set() for _ in rows]; pc_exons=[set() for _ in rows]; pct_exons=[set() for _ in rows]
    genes={}; gtf_counts=collections.Counter(); retained_intervals=collections.Counter(); metadata=[]
    gene_re=re.compile(r'(?:^|; )gene_id "([^"]+)"')
    type_re=re.compile(r'(?:^|; )gene_type "([^"]+)"')
    name_re=re.compile(r'(?:^|; )gene_name "([^"]+)"')
    tx_re=re.compile(r'(?:^|; )transcript_type "([^"]+)"')
    h=hashlib.sha256(); known_boundaries=[]
    with GTF.open('rb') as f:
        for raw in f:
            h.update(raw)
            if raw.startswith(b'#'):
                metadata.append(raw.decode().rstrip());continue
            q=raw.decode().rstrip('\n').split('\t',8)
            feature=q[2];gtf_counts[feature]+=1
            if feature not in ('gene','exon'): continue
            chrom,start,end=q[0],int(q[3]),int(q[4]); attrs=q[8]
            gid=gene_re.search(attrs).group(1); gtype=type_re.search(attrs).group(1)
            if feature=='gene':
                name=name_re.search(attrs).group(1)
                assert gid not in genes, f'Duplicate versioned gene ID {gid}'
                genes[gid]={'gene_id':gid,'gene_symbol':name,'gene_type':gtype,'chrom':chrom,'start1':start,'end1':end,'strand':q[6]}
                if name in ('BAP1','RASA1','LATS1','LATS2','APOE','ACTB','HLA-B'):
                    ps=[start-1,start,end,end+1]
                    observed=[ps[i] for i in indices_in_closed_interval(ps,start,end)]
                    assert observed==[start,end]
                    known_boundaries.append({'gene_id':gid,'symbol':name,'chrom':chrom,'start1':start,'end1':end,'test_positions':ps,'included':observed,'expected':[start,end],'status':'PASS'})
                destination=spans
            else:
                if gtype!='protein_coding': continue
                destination=pc_exons
            if chrom not in positions: continue
            indices=indices_in_closed_interval(positions[chrom],start,end)
            if not indices: continue
            retained_intervals[feature]+=1
            tx_pc=feature=='exon' and tx_re.search(attrs).group(1)=='protein_coding'
            for j in indices:
                i=rowids[chrom][j];destination[i].add(gid)
                if tx_pc:pct_exons[i].add(gid)
    assert len(known_boundaries)>=7
    checks['known_GENCODE_gene_boundaries']=known_boundaries
    # Exons must sit inside a corresponding actual gene span; preserve all overlapping genes.
    for i in range(len(rows)):
        assert pc_exons[i] <= spans[i] and pct_exons[i] <= pc_exons[i]
        assert all(genes[g]['gene_type']=='protein_coding' for g in pc_exons[i])
    checks['exon_gene_subset_checks']='PASS for all 165782 loci'
    category=collections.Counter();bychr=collections.defaultdict(collections.Counter);byscope=collections.defaultdict(collections.Counter)
    gene_counts=collections.Counter(); unique_gene_counts=collections.Counter(); pc_gene_counts=collections.Counter()
    fields=list(rows[0])+['gene_overlap_category','n_overlapping_genes','overlap_gene_ids','overlap_gene_symbols','overlap_gene_types','unique_gene_id','unique_gene_symbol','unique_gene_type','protein_coding_gene_ids','protein_coding_gene_symbols','pc_gene_exon_gene_ids','pc_gene_exon_gene_symbols','pc_transcript_exon_gene_ids','pc_transcript_exon_gene_symbols']
    target=OUT/'public-panel-gene-annotations.tsv'
    def joined(ids,field):return ';'.join(genes[g][field] for g in ids) or '.'
    with target.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for i,r in enumerate(rows):
            gids=sorted(spans[i]);pcids=[g for g in gids if genes[g]['gene_type']=='protein_coding']; ex=sorted(pc_exons[i]);txex=sorted(pct_exons[i])
            cat='intergenic' if not gids else 'unique_gene' if len(gids)==1 else 'multi_gene'
            cats=[cat,'all_loci']+(['in_pc_gene_exon'] if ex else [])+(['in_pc_transcript_exon'] if txex else [])+(['overlaps_protein_coding_gene'] if pcids else [])
            if len(pcids)==1:cats.append('exactly_one_pc_gene_regardless_other_gene_overlap')
            for key in cats:category[key]+=1;bychr[r['chrom']][key]+=1;byscope[r['scope']][key]+=1
            gene_counts.update(gids);pc_gene_counts.update(pcids)
            if len(gids)==1:unique_gene_counts[gids[0]]+=1
            out=dict(r,gene_overlap_category=cat,n_overlapping_genes=len(gids),overlap_gene_ids=joined(gids,'gene_id'),overlap_gene_symbols=joined(gids,'gene_symbol'),overlap_gene_types=joined(gids,'gene_type'),unique_gene_id=gids[0] if len(gids)==1 else '.',unique_gene_symbol=genes[gids[0]]['gene_symbol'] if len(gids)==1 else '.',unique_gene_type=genes[gids[0]]['gene_type'] if len(gids)==1 else '.',protein_coding_gene_ids=joined(pcids,'gene_id'),protein_coding_gene_symbols=joined(pcids,'gene_symbol'),pc_gene_exon_gene_ids=joined(ex,'gene_id'),pc_gene_exon_gene_symbols=joined(ex,'gene_symbol'),pc_transcript_exon_gene_ids=joined(txex,'gene_id'),pc_transcript_exon_gene_symbols=joined(txex,'gene_symbol'))
            w.writerow(out)
    with (OUT/'overlapping-gene-inventory.tsv').open('w',newline='') as f:
        fields2=['gene_id','gene_symbol','gene_type','chrom','start1','end1','strand','panel_overlap_loci','panel_uniquely_assigned_loci']
        w=csv.DictWriter(f,fieldnames=fields2,delimiter='\t',lineterminator='\n');w.writeheader()
        for gid in sorted(gene_counts):w.writerow(dict(genes[gid],panel_overlap_loci=gene_counts[gid],panel_uniquely_assigned_loci=unique_gene_counts[gid]))
    summary={'total_loci':len(rows),'categories':dict(category),'distinct_overlapping_gene_ids':len(gene_counts),'distinct_gene_ids_with_at_least_one_unique_gene_locus':len(unique_gene_counts),'distinct_overlapping_protein_coding_gene_ids':len(pc_gene_counts),'distinct_protein_coding_gene_ids_with_at_least_one_unique_gene_locus':sum(genes[g]['gene_type']=='protein_coding' for g in unique_gene_counts),'by_chromosome':{k:dict(v) for k,v in sorted(bychr.items())},'by_scope':{k:dict(v) for k,v in sorted(byscope.items())},'annotation_gene_count_all_contigs':len(genes),'GTF_feature_counts':dict(gtf_counts),'GTF_intervals_overlapping_panel':dict(retained_intervals),'scope':'Public panel only. Counts describe locus annotations, not patient coverage, observed variants, genotypes, expression or sample concordance.','multi_gene_policy':'Any two or more versioned gene IDs with overlapping gene spans are multi_gene even if only one is protein_coding. No unique assignment is imputed.','exon_policy':'pc_gene_exon includes exons from any transcript of a protein_coding gene; pc_transcript_exon restricts transcript_type to protein_coding. Both include UTR exonic bases and do not mean CDS membership.','intergenic_definition':'No overlapping GENCODE37 gene feature on this exact primary-reference chromosome; not a claim of absent regulatory function.'}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    (OUT/'boundary-validation.json').write_text(json.dumps(checks,indent=2))
    provenance={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script':str(Path(__file__).resolve()),'script_sha256':sha(Path(__file__)),'input_panel':{'path':str(PANEL),'sha256':sha(PANEL),'bytes':PANEL.stat().st_size,'rows':len(rows),'order':'Output preserves frozen input row order'},'annotation':{'path':str(GTF),'sha256_uncompressed':h.hexdigest(),'bytes':GTF.stat().st_size,'source_manifest_path':str(SOURCE),'source_manifest_sha256':sha(SOURCE),'original_download':json.loads(SOURCE.read_text()),'metadata':metadata,'build':'GRCh38 primary assembly, GENCODE37 / Ensembl103'},'coordinate_conventions':'Input pos1 and GTF start/end are 1-based inclusive. A locus overlaps iff start1 <= pos1 <= end1 on the same exact contig name. No nearest-gene imputation, no strand filter, no extension/flanking window.','patient_data_read':False,'outputs':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in [target,OUT/'overlapping-gene-inventory.tsv',OUT/'summary.json',OUT/'boundary-validation.json']},'elapsed_seconds':round(time.time()-started,3),'peak_RSS_bytes_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2))
    print(json.dumps({'categories':summary['categories'],'distinct_unique_assignment_genes':summary['distinct_gene_ids_with_at_least_one_unique_gene_locus'],'elapsed_seconds':provenance['elapsed_seconds'],'RSS_bytes':provenance['peak_RSS_bytes_macOS'],'outputs':provenance['outputs']},indent=2))

if __name__=='__main__':main()

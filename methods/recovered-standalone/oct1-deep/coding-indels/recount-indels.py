"""Fixed-window reference/alternate haplotype recount; exploratory only."""
from pathlib import Path
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
import json,urllib.request,time,csv
import pysam
ROOT=Path(__file__).resolve().parents[3]; OUT=Path(__file__).resolve().parent
SRC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
REFDIR=OUT/'reference-windows'; REFDIR.mkdir(exist_ok=True)
CAND=json.loads((OUT/'strict-frameshift-shortlist.json').read_text())['candidates']

def reference(v):
    p=v['pos1']-1;start=max(0,p-80);end=p+len(v['ref'])+80
    fn=REFDIR/f'{v["chrom"]}-{start}-{end}.json'
    if fn.exists(): return json.loads(fn.read_text())
    url=f'https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={v["chrom"]};start={start};end={end}'
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url,timeout=40) as r: data=json.load(r)
            assert len(data['dna'])==end-start
            data['source_url']=url;fn.write_text(json.dumps(data)+'\n');return data
        except Exception:
            if attempt==2:raise
            time.sleep(1+attempt)

def prepare(v):
    ref=reference(v);seq=ref['dna'].upper(); p=v['pos1']-1
    i=p-ref['start'];assert seq[i:i+len(v['ref'])]==v['ref'],v
    start=p-12;end=p+len(v['ref'])+12
    rseq=seq[start-ref['start']:end-ref['start']]
    aseq=rseq[:12]+v['alt']+rseq[12+len(v['ref']):]
    return start,end,rseq,aseq,ref

def read_haplotype(read,start,end):
    q=0;r=read.reference_start;pieces=[];qindexes=[];coverleft=False;coverright=False
    for op,n in read.cigartuples:
        if op in (0,7,8):
            a=max(start,r);b=min(end,r+n)
            if a<b:
                qa=q+a-r;qb=q+b-r;pieces.append(read.query_sequence[qa:qb]);qindexes.extend(range(qa,qb))
                coverleft|=a==start;coverright|=b==end
            q+=n;r+=n
        elif op==1:
            if start<r<end:pieces.append(read.query_sequence[q:q+n]);qindexes.extend(range(q,q+n))
            q+=n
        elif op==2:r+=n
        elif op==3:
            if r<end and r+n>start:return None,'spliced_across_window'
            r+=n
        elif op==4:q+=n
    if not coverleft or not coverright:return None,'incomplete_window'
    if not qindexes or read.query_qualities is None or min(read.query_qualities[i] for i in qindexes)<20:return None,'low_base_quality'
    return ''.join(pieces),'assessable'

def recount(bam,v,setup):
    start,end,ref,alt,_=setup; counts=Counter();calls=defaultdict(set);strands=defaultdict(Counter);clean=defaultdict(set);ends=defaultdict(set)
    for read in bam.fetch(v['chrom'],start,end):
        if read.flag&(4|256|512|1024|2048) or read.mapping_quality<20:counts['excluded_flag_or_mapq']+=1;continue
        h,reason=read_haplotype(read,start,end)
        if h is None:counts[reason]+=1;continue
        call='alt' if h==alt else ('ref' if h==ref else 'other_haplotype')
        counts[call]+=1;name=(read.get_tag('RG') if read.has_tag('RG') else '',read.query_name);calls[name].add(call)
        strands[call]['reverse' if read.is_reverse else 'forward']+=1
        ends[call].add((read.reference_start,read.reference_end,read.is_reverse))
        if not any(op==4 for op,n in read.cigartuples):clean[call].add(name)
    fragments=Counter(next(iter(vals)) if len(vals)==1 else 'discordant' for vals in calls.values())
    den=fragments['alt']+fragments['ref']
    return {'reads':dict(counts),'paired_names':dict(fragments),'alt_fraction_ref_alt_only':fragments['alt']/den if den else None,'orientations':{k:dict(v) for k,v in strands.items()},'no_softclip_names':{k:len(v) for k,v in clean.items()},'distinct_endpoints':{k:len(v) for k,v in ends.items()}}

def main():
    started=time.time();setups=list(ThreadPoolExecutor(max_workers=4).map(prepare,CAND));result=[]
    handles={kind:pysam.AlignmentFile(str(SRC/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai')) for kind in ['DNA','RNA']}
    for i,(v,setup) in enumerate(zip(CAND,setups)):
        row=dict(v);row['reference_source']=setup[4]['source_url'];row['reference_start0']=setup[0];row['reference_end0']=setup[1];row['reference_haplotype']=setup[2];row['alternate_haplotype']=setup[3]
        row['DNA']=recount(handles['DNA'],v,setup)
        dna=row['DNA']['paired_names'];fraction=row['DNA']['alt_fraction_ref_alt_only']
        # RNA is recounted even for DNA-low-fraction candidates so negative
        # filtering is transparent and positive/negative controls remain visible.
        row['RNA']=recount(handles['RNA'],v,setup)
        row['DNA_threshold_pass']=dna.get('alt',0)>=10 and fraction is not None and fraction>=.05 and min(row['DNA']['orientations'].get('alt',{}).get(k,0) for k in ['forward','reverse'])>=3
        row['RNA_support_ge3_names']=row['RNA']['paired_names'].get('alt',0)>=3
        result.append(row)
        print(i+1,','.join(v['genes']),v['pos1'],v['ref'],v['alt'],'DNA',dna,'RNA',row['RNA']['paired_names'],flush=True)
    for b in handles.values():b.close()
    obj={'method':'Fixed 12-reference-base flanks around each raw CIGAR VCF-style allele. Full contiguous local haplotype must exactly match public hg38 reference or edited alternate, with every included baseQ>=20 and MAPQ>=20, primary/nonduplicate/QCpassing. Paired names collapsed and conflicting mates excluded. Repeats/nearby variants can reduce assessability; no somatic classification or clinical sensitivity claim. Fractions use only exact ref+alt names. RNA no duplicate flags: names are not independent UMI molecules.','candidates':result,'elapsed_seconds':time.time()-started}
    (OUT/'coding-indel-haplotype-recount.json').write_text(json.dumps(obj,indent=2)+'\n')
    fields=['gene','chrom','pos1','ref','alt','DNA_alt','DNA_ref','DNA_fraction','RNA_alt','RNA_ref','DNA_threshold_pass','RNA_support_ge3_names']
    with (OUT/'coding-indel-haplotype-recount.tsv').open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader()
        for r in result:
            w.writerow({'gene':';'.join(r['genes']),'chrom':r['chrom'],'pos1':r['pos1'],'ref':r['ref'],'alt':r['alt'],'DNA_alt':r['DNA']['paired_names'].get('alt',0),'DNA_ref':r['DNA']['paired_names'].get('ref',0),'DNA_fraction':r['DNA']['alt_fraction_ref_alt_only'],'RNA_alt':r['RNA']['paired_names'].get('alt',0),'RNA_ref':r['RNA']['paired_names'].get('ref',0),'DNA_threshold_pass':r['DNA_threshold_pass'],'RNA_support_ge3_names':r['RNA_support_ge3_names']})
    print('Completed',len(result),'DNA/RNA candidates; seconds',time.time()-started,flush=True)

if __name__=='__main__':main()

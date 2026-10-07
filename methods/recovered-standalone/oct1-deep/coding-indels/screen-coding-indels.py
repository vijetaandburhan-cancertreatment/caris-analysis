"""Bounded coding-indel discovery from supplied tumor BAM; not a clinical caller.

No independent somatic classification, calibrated allele fractions, or negative
clinical claims. A downstream local-reference/haplotype recount is required.
"""
from pathlib import Path
from collections import defaultdict, Counter
import argparse, gzip, json, re, time, resource
import pysam

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SRC = Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
GTF = ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz'
INDEX = ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai'
BIN = 10000

def annotations():
    saved = OUT/'selected-coding-transcripts.json'
    if saved.exists():
        return json.loads(saved.read_text())
    tx = {}
    with gzip.open(GTF, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            a = line.rstrip('\n').split('\t')
            if len(a) != 9 or a[2] != 'CDS': continue
            fields = defaultdict(list)
            for k,v in re.findall(r'(\w+) "([^"]+)"',a[8]): fields[k].append(v)
            if fields['gene_type'] != ['protein_coding']: continue
            ident = fields['transcript_id'][0]
            if ident not in tx:
                tags=fields['tag']
                priority = 0 if 'MANE_Select' in tags else (1 if 'appris_principal_1' in tags else (2 if any(t.startswith('appris_principal') for t in tags) else 3))
                tx[ident]={'transcript':ident,'gene':fields['gene_name'][0],'gene_id':fields['gene_id'][0], 'chrom':a[0], 'strand':a[6], 'tags':tags, 'priority':priority,'CDS':[]}
            tx[ident]['CDS'].append([int(a[3])-1,int(a[4])])
    by_gene=defaultdict(list)
    for d in tx.values():
        d['CDS']=sorted(d['CDS']); d['coding_bases']=sum(e-s for s,e in d['CDS'])
        by_gene[d['gene_id']].append(d)
    result=[]
    for gene,candidates in by_gene.items():
        result.append(sorted(candidates,key=lambda d:(d['priority'],-d['coding_bases'],d['transcript']))[0])
    saved.write_text(json.dumps(result,indent=2)+'\n')
    return result

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--region'); ap.add_argument('--prefix',default='all'); args=ap.parse_args()
    started=time.time(); ann=annotations(); bins=defaultdict(list)
    for i,g in enumerate(ann):
        for s,e in g['CDS']:
            for b in range(s//BIN,(e-1)//BIN+1): bins[(g['chrom'],b)].append((s,e,i))
    stats=Counter(); events={}; lasttime=time.time()
    def genes_at(chrom,s,e):
        found=set()
        for b in range(max(0,s)//BIN,max(0,e-1)//BIN+1):
            for x,y,i in bins.get((chrom,b),[]):
                if x<e and y>s: found.add(i)
        return found
    bam=pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(INDEX),threads=2)
    stream=bam.fetch(region=args.region) if args.region else bam.fetch(until_eof=True)
    for r in stream:
        stats['records_seen']+=1
        if stats['records_seen']%2000000==0:
            report={'records_seen':stats['records_seen'],'retained_event_keys':len(events),'elapsed_seconds':round(time.time()-started),'max_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'last_coordinate':[r.reference_name,r.reference_start]}
            (OUT/(args.prefix+'-progress.json')).write_text(json.dumps(report)+'\n')
            print(report,flush=True)
        if r.flag & (4|256|512|1024|2048) or r.mapping_quality<20: continue
        cig=r.cigartuples
        if not cig or not any(op in (1,2) and n<=50 for op,n in cig): continue
        stats['qualifying_reads_with_indel']+=1
        if any(op==3 for op,n in cig):
            stats['skipped_DNA_reference_skip_reads']+=1; continue
        if not r.has_tag('MD'):
            stats['missing_MD_indel_reads']+=1; continue
        quals=r.query_qualities; seq=r.query_sequence
        if quals is None or seq is None: continue
        chrom=r.reference_name; refpos=r.reference_start; qpos=0
        refseq=None; no_softclip=not any(op==4 for op,n in cig)
        for j,(op,n) in enumerate(cig):
            if op in (0,7,8): refpos+=n; qpos+=n; continue
            if op==4: qpos+=n; continue
            if op not in (1,2): continue
            alen=n if op==1 else 0
            gene_ids=genes_at(chrom,refpos-1,refpos+(n if op==2 else 1)) if 1<=n<=50 else set()
            if gene_ids and j>0 and j+1<len(cig) and cig[j-1][0] in (0,7,8) and cig[j+1][0] in (0,7,8) and cig[j-1][1]>=8 and cig[j+1][1]>=8:
                flankquals=quals[qpos-8:qpos+alen+8]
                if len(flankquals)==16+alen and min(flankquals)>=20:
                    if refseq is None:
                        try: refseq=r.get_reference_sequence().upper()
                        except Exception: stats['MD_reconstruction_errors']+=1; refseq=''
                    offset=refpos-r.reference_start
                    if refseq and 0<offset<len(refseq):
                        ref=refseq[offset-1:offset+(n if op==2 else 0)]
                        alt=refseq[offset-1]+seq[qpos:qpos+n] if op==1 else refseq[offset-1]
                        if set(ref+alt)<=set('ACGT') and len(ref)==(n+1 if op==2 else 1):
                            key=(chrom,refpos,ref,alt)
                            if key not in events: events[key]={'reads':0,'fragments':{},'endpoints':set(),'genes':set(),'mapq60_reads':0,'clean_reads':0}
                            d=events[key]; d['reads']+=1; d['genes'].update(gene_ids)
                            name=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name)
                            d['fragments'][name]=d['fragments'].get(name,0)|(2 if r.is_reverse else 1)
                            if len(d['endpoints'])<200: d['endpoints'].add((r.reference_start,r.reference_end,r.is_reverse))
                            d['mapq60_reads']+=r.mapping_quality>=60;d['clean_reads']+=no_softclip
                            stats['quality_filtered_coding_indel_observations']+=1
            if op==1:qpos+=n
            else:refpos+=n
    bam.close()
    rows=[]
    for (chrom,pos,ref,alt),d in events.items():
        n=len(d['fragments'])
        if n<5:continue
        row={'chrom':chrom,'pos1':pos,'ref':ref,'alt':alt,'net_length_change':len(alt)-len(ref),'frameshift_length':(len(alt)-len(ref))%3!=0,'alternate_reads':d['reads'],'alternate_paired_names':n,'forward_supported_names':sum(v&1!=0 for v in d['fragments'].values()),'reverse_supported_names':sum(v&2!=0 for v in d['fragments'].values()),'distinct_alignment_endpoints_capped200':len(d['endpoints']),'mapq60_reads':d['mapq60_reads'],'no_softclip_reads':d['clean_reads'],'genes':[ann[i]['gene'] for i in sorted(d['genes'])],'selected_transcripts':[ann[i]['transcript'] for i in sorted(d['genes'])]}
        rows.append(row)
    rows.sort(key=lambda r:(r['frameshift_length'],r['alternate_paired_names']),reverse=True)
    result={'method':'Research discovery only: one selected GENCODE37 protein-coding transcript per gene (MANE then APPRIS then longest), coding overlap, primary MAPQ>=20, exclude duplicate/QCfail/secondary/supplementary/unmapped, CIGAR indels 1..50 bp, immediately adjacent >=8 aligned bases on each side with baseQ>=20, inserted baseQ>=20, reference reconstructed from BAM MD. Report >=5 alternate paired names. Raw CIGAR alleles not normalized; no reference depth denominator or somatic classification. Follow-up reference/haplotype audit required.','stats':dict(stats),'selected_genes':len(ann),'event_keys_with_any_support':len(events),'reported_candidates':len(rows),'elapsed_seconds':time.time()-started,'max_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'candidates':rows}
    (OUT/(args.prefix+'-coding-indel-discovery.json')).write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ('method','candidates')},flush=True)

if __name__=='__main__':main()

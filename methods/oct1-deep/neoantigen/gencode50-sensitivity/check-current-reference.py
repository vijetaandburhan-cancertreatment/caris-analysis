"""Bounded public GENCODE download and exhaustive local peptide window scan.

No patient sequence is sent over the network. Searches every equal-length
window, with no seeds, no gaps and no undocumented affinity inference.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
import csv, gzip, hashlib, json, resource, time, urllib.request
import numpy as np
from Bio import SeqIO

OUT = Path(__file__).resolve().parent
PARENT = OUT.parent
URL = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.pc_translations.fa.gz'
CURRENT = OUT / 'gencode.v50.pc_translations.fa.gz'
CAP = 30_000_000
QUERIES = {'IEERKGLYL':2, 'EERKGLYL':2, 'RRFSSLQT':0}

def digest(path, algo='sha256'):
    h=hashlib.new(algo)
    with path.open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()

def get_public_reference():
    started=datetime.now(timezone.utc).isoformat()
    head=urllib.request.Request(URL,method='HEAD')
    with urllib.request.urlopen(head,timeout=45) as r:
        headers=dict(r.headers);size=int(r.headers.get('Content-Length','0'))
    if size <=0 or size>CAP:raise RuntimeError(f'Reference size {size} exceeds permitted scope or unavailable')
    if not CURRENT.exists():
        tmp=CURRENT.with_suffix(CURRENT.suffix+'.partial')
        total=0
        with urllib.request.urlopen(URL,timeout=60) as r,tmp.open('wb') as f:
            while block:=r.read(1024*1024):
                total+=len(block)
                if total>CAP:raise RuntimeError('Download exceeds30MB cap')
                f.write(block)
        if total!=size:raise RuntimeError(f'Incomplete public download {total} != {size}')
        tmp.rename(CURRENT)
    if CURRENT.stat().st_size!=size:raise RuntimeError('Cached reference size differs from current server header')
    manifest={'download_started_utc':started,'release':'GENCODE human50','release_month_from_official_history':'2026-06','assembly':'GRCh38.p14','official_page':'https://www.gencodegenes.org/human/','history_page':'https://www.gencodegenes.org/human/releases.html','url':URL,'compressed_size_bytes':size,'sha256':digest(CURRENT),'md5':digest(CURRENT,'md5'),'server_headers':headers,'scope':'Official protein-coding transcript translations, all regions including scaffolds/patches/alternate loci; biotypes include protein_coding, nonsense_mediated_decay, non_stop_decay, IG/TR, polymorphic_pseudogene, protein_coding_LoF. This is a reference translation collection, not an empirical normal-tissue proteomics atlas.'}
    (OUT/'public-reference-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest

def exhaustive_scan(path,label):
    started=time.time();hits=defaultdict(list);records=0;residues=0;windows={q:0 for q in QUERIES};query_arrays={q:np.frombuffer(q.encode('ascii'),dtype=np.uint8) for q in QUERIES}
    with gzip.open(path,'rt') as f:
        for record in SeqIO.parse(f,'fasta'):
            records+=1;seq=str(record.seq);residues+=len(seq)
            arr=np.frombuffer(seq.encode('ascii'),dtype=np.uint8)
            for q,threshold in QUERIES.items():
                n=len(q);count=len(seq)-n+1
                if count<=0:continue
                windows[q]+=count
                mismatches=np.zeros(count,dtype=np.uint8)
                for offset,aa in enumerate(query_arrays[q]):mismatches+=arr[offset:offset+count]!=aa
                for s in np.flatnonzero(mismatches<=threshold):
                    s=int(s);p=seq[s:s+n]
                    hits[(q,p)].append({'reference_header':record.description,'start_aa1':s+1,'end_aa1':s+n,'substitution_positions1':[i+1 for i,(a,b) in enumerate(zip(q,p)) if a!=b],'left_flank5':seq[max(0,s-5):s],'right_flank5':seq[s+n:s+n+5]})
            if records%50000==0:print(label,records,'records',round(time.time()-started,1),'seconds',flush=True)
    rows=[{'query':q,'reference_peptide':p,'substitutions':sum(a!=b for a,b in zip(q,p)),'reference_occurrences':len(where),'occurrences':where} for (q,p),where in sorted(hits.items())]
    summary={q:{'maximum_hamming_searched':QUERIES[q],'windows_examined':windows[q],'distinct_peptides_by_substitutions':{str(n):sum(r['query']==q and r['substitutions']==n for r in rows) for n in range(QUERIES[q]+1)},'reference_occurrences_by_substitutions':{str(n):sum(r['reference_occurrences'] for r in rows if r['query']==q and r['substitutions']==n) for n in range(QUERIES[q]+1)}} for q in QUERIES}
    result={'label':label,'compressed_path':str(path),'compressed_sha256':digest(path),'records':records,'residues':residues,'elapsed_seconds':round(time.time()-started,3),'summary':summary,'results':rows,'method':'Exhaustive unseeded equal-length window scan with NumPy uint8 mismatch counts; stream one FASTA record at a time. Hamming<=2 for BAP1 9/8mer, exact only for optional RASA1 8mer. No gap/length-tolerant matching, population-variant proteome, noncanonical translation, PTM or TCR cross-reactivity modeling.'}
    (OUT/f'{label}-window-results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def main():
    started=time.time();manifest=get_public_reference()
    old=exhaustive_scan(PARENT/'gencode.v37.pc_translations.fa.gz','gencode37')
    new=exhaustive_scan(CURRENT,'gencode50')
    prior=json.loads((PARENT/'bap1-normal-sequence-neighbors.json').read_text())
    oldpairs={(r['query'],r['reference_peptide']) for r in old['results'] if r['query']!='RRFSSLQT'}
    seededpairs={(r['query'],r['normal_reference_peptide']) for r in prior['results']}
    assert oldpairs==seededpairs,'Independent exhaustive G37 scan disagrees with previous seeded search'
    comparisons={}
    for q in QUERIES:
        a={r['reference_peptide'] for r in old['results'] if r['query']==q}
        b={r['reference_peptide'] for r in new['results'] if r['query']==q}
        comparisons[q]={'added_distinct_peptides':sorted(b-a),'removed_distinct_peptides':sorted(a-b),'shared_distinct_count':len(a&b),'gencode37':old['summary'][q],'gencode50':new['summary'][q]}
    report={'finished_utc':datetime.now(timezone.utc).isoformat(),'public_reference':manifest,'records_scanned':{'GENCODE37':old['records'],'GENCODE50':new['records']},'G37_independent_unseeded_agrees_with_existing_seeded_method':True,'comparison':comparisons,'peak_RSS_bytes_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'elapsed_seconds':round(time.time()-started,2),'interpretation_limit':'No exact reference peptide match is not proof of tumor specificity, translation, HLA display, safety or T-cell recognition. Normal population variants and unannotated/noncanonical ORFs remain outside this reference. Similar sequences are possible experimental controls, not validated off-targets.'}
    (OUT/'findings.json').write_text(json.dumps(report,indent=2)+'\n')
    with (OUT/'gencode50-neighbors.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['query','reference_peptide','substitutions','reference_occurrences','example_header'])
        for r in new['results']:w.writerow([r['query'],r['reference_peptide'],r['substitutions'],r['reference_occurrences'],r['occurrences'][0]['reference_header']])
    print(json.dumps({k:v for k,v in report.items() if k!='public_reference'},indent=2),flush=True)

if __name__=='__main__':main()

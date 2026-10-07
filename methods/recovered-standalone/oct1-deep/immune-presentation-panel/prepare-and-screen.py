"""Bounded coding screen of antigen-presentation/interferon machinery.

No inference of intact pathway function or immunotherapy efficacy from absence
of a candidate. HLA genes deliberately excluded from conventional BAM calls.
"""
from pathlib import Path
import json,gzip,re,runpy
B=Path(__file__).resolve().parent;ROOT=B.parents[2]
GENES=['B2M','TAP1','TAP2','TAPBP','JAK1','JAK2','IFNGR1','IFNGR2','STAT1','IRF1','NLRC5','PSMB8','PSMB9','CIITA','RFX5','RFXANK','RFXAP']
selected={x['gene']:x for x in json.loads((ROOT/'work/oct1-deep/coding-indels/selected-coding-transcripts.json').read_text()) if x['gene'] in GENES}
genes={}
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
    for line in f:
        if line.startswith('#'):continue
        a=line.rstrip().split('\t')
        if a[2]!='gene' or a[0] not in ['chr'+str(i) for i in range(1,23)]:continue
        m=re.search(r'gene_name "([^"]+)"',a[8])
        if m and m[1] in selected:
            g=m[1];s=selected[g]
            genes[g]={'gene':g,'chrom':a[0],'strand':a[6],'gene_start0':int(a[3])-1,'gene_end0':int(a[4]),'selected_transcript':s['transcript'],'selected_CDS':sorted(s['CDS'])}
assert set(genes)==set(GENES)
(B/'target-annotation.json').write_text(json.dumps({'genes':genes},indent=2)+'\n')
source=(ROOT/'work/oct1-deep/extended-panel/screen.py').read_text()
source=source.replace('Bounded candidate discovery in 13 named genes; no clinical calling.','Bounded candidate discovery in 17 antigen-presentation/interferon genes; no clinical calling.')
old="GENES=['NF2','SMARCA4','SMARCB1','SETD2','TP53','BRAF','KRAS','NRAS','EGFR','ERBB2','RB1','PTEN','PBRM1']"
assert old in source
source=source.replace(old,'GENES='+repr(GENES))
source=source.replace('Strict read base-depth screen across MANE CDS+/-4bp','Strict read base-depth screen across selected CDS+/-4bp (16 MANE transcripts; NLRC5 APPRIS principal_4)')
(B/'screen.py').write_text(source)
runpy.run_path(str(B/'screen.py'),run_name='__main__')

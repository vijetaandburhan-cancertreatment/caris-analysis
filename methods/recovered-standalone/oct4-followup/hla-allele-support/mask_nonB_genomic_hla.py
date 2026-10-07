from pathlib import Path
import sys,gzip,json,re,time,resource
from collections import defaultdict,Counter
from Bio import SeqIO
ROOT=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-hla-allele-support');sys.path.insert(0,str(ROOT/'deps'));import ahocorasick
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with gzip.open(ROOT/'common-haplotype-markers.json.gz','rt') as f:j=json.load(f)
ms=[m for m in j['markers'] if m['B_gene_exclusive'] and m['passes_transcript_paralog_mask']];idx=defaultdict(set)
for i,m in enumerate(ms):
 for s in {m['sequence'],rc(m['sequence'])}:idx[s].add(i)
a=ahocorasick.Automaton()
for s,ids in idx.items():a.add_word(s,tuple(ids))
a.make_automaton();hits=defaultdict(set);ct=Counter();t0=time.time()
for r in SeqIO.parse('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools/arcasHLA/dat/IMGTHLA/hla.dat','imgt'):
 m=re.search(r'^(?:HLA-)?([A-Za-z0-9]+)\*',r.description)
 if not m:ct['unparsed_gene_records']+=1;continue
 g=m.group(1)
 if g=='B':ct['B_records_not_paralog_screen']+=1;continue
 ct['nonB_genomic_records']+=1;ct['nonB_genomic_bases']+=len(r.seq)
 for i in {i for p,ids in a.iter(str(r.seq).upper()) for i in ids}:hits[i].add(g)
out={'status':'passed','counts':dict(ct),'new_nonB_matches':[{'marker_index_filtered':i,'allele':ms[i]['allele'],'k':ms[i]['k'],'sequence':ms[i]['sequence'],'genes':sorted(v)} for i,v in hits.items()],'newly_masked_markers':len(hits),'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'purpose':'Screen every non-HLA-B genomic record, including records without CDS, for paralog/pseudogene sequence matches. Complements all-CDS and whole-transcriptome marker masks.'}
(ROOT/'nonB-genomic-mask.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='new_nonB_matches'}),flush=True)

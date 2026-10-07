"""Reconstruct descriptive read-context metrics from the preserved indexed extraction."""
from pathlib import Path
import json,gzip,collections
P=Path(__file__).parent
rows=[json.loads(x)for x in gzip.open(P/'all-site-read-records.jsonl.gz','rt')];out={}
for site,base in {'chr2:189043211':'A','chr3:195568834':'C','chr5:443259':'A','chr5:80654917':'G','chr8:100705591':'G','chr8:100705604':'C'}.items():
 rr=[r for r in rows if r['site']==site and r['assay']=='RNA'and r['baseline']];ns=collections.defaultdict(set)
 for r in rr:ns[(r['name'],r['RG'])].add(r['observation']['base'])
 rr=[r for r in rr if r['observation']['base']==base and ns[(r['name'],r['RG'])]=={base}]
 out[site]={'unexpected':base,'reads':len(rr),'all_qualifying_name_units_including_conflicts':len(ns),'unexpected_names':len({r['name']for r in rr}),'aligned_left_hist':dict(collections.Counter(r['observation']['aligned_query_left']for r in rr)),'within5_aligned_edge_reads':sum(min(r['observation']['aligned_query_left'],r['observation']['aligned_query_right'])<5 for r in rr),'softclip_reads':sum(any(e['op']=='S'for e in r['observation']['near_CIGAR_events'])for r in rr),'nearby_events':dict(collections.Counter((e['op']+':'+str(e['length'])+':'+str(e['rpos0'])+':'+str(e.get('rend0','')))for r in rr for e in r['observation']['near_CIGAR_events']if e['distance_genomic']<50))}
g=collections.defaultdict(lambda:collections.defaultdict(set));same=collections.defaultdict(dict)
for r in rows:
 if r['assay']=='RNA'and r['site']in('chr8:100705591','chr8:100705604')and r['baseline']:
  g[(r['name'],r['RG'])][r['site']].add(r['observation']['base']);k=(r['name'],r['RG'],r['flag'],r['start0'],r['CIGAR'],r['sequence']);same[k][r['site']]=r['observation']['base']
ct=collections.Counter();sam=collections.defaultdict(set)
for name,v in g.items():
 if len(v)==2 and all(len(s)==1 for s in v.values()):ct[next(iter(v['chr8:100705591']))+next(iter(v['chr8:100705604']))]+=1
for k,v in same.items():
 if len(v)==2:sam[v['chr8:100705591']+v['chr8:100705604']].add(k[:2])
out['PABPC1_phase']={'unambiguous_names_covering_both':dict(ct),'names_with_one_same_alignment_covering_both':{k:len(v)for k,v in sam.items()},'note':'Same-query and same-record observations, not UMI-independent molecules or genomic phase proof.'}
# Compare reconstructed values to the already generated output without overwriting it.
old=json.loads((P/'specific-context-metrics.json').read_text());new=json.loads(json.dumps(out));assert old==new
print('Preserved metrics independently reproducible from extracted records')

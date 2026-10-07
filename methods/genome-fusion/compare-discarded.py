"""Read-only exact row-multiset comparison and explicit characterization of bounded differences."""
import pathlib,collections,hashlib,json,time
P=pathlib.Path(__file__).resolve().parent;S=P/'patient-D8/fusions.discarded.tsv';D=P/'diagnostics/patient-D8-counter-rerun/fusions.discarded.tsv';counts=[];t=time.time()
for f in [S,D]:
 c=collections.Counter();n=0
 with f.open('rb') as h:
  header=next(h)
  for row in h:c[hashlib.sha256(row).digest()]+=1;n+=1
 counts.append((c,n,header))
a,an,ah=counts[0];b,bn,bh=counts[1];left=a-b;right=b-a
out={'stock_rows':an,'diagnostic_rows':bn,'headers_identical':ah==bh,'exact_row_multisets_identical':a==b,'stock_only_row_instances':sum(left.values()),'diagnostic_only_row_instances':sum(right.values()),'stock_only_unique_rows':len(left),'diagnostic_only_unique_rows':len(right),'elapsed_seconds':time.time()-t}
for f,key,want in [(S,'stock_examples',set(left)),(D,'diagnostic_examples',set(right))]:
 out[key]=[]
 with f.open('rb') as h:
  cols=next(h).decode().rstrip('\n').lstrip('#').split('\t')
  for row in h:
   if hashlib.sha256(row).digest() in want:
    out[key].append(dict(zip(cols,row.decode().rstrip('\n').split('\t'))))
    if len(out[key])>=10:break
differences=[]
for s in out['stock_examples']:
 matches=[r for r in out['diagnostic_examples'] if all(s[k]==r[k] for k in s if k!='gene1')]
 if len(matches)==1:
  r=matches[0];differences.append({'breakpoints':[s['breakpoint1'],s['breakpoint2']],'gene2':s['gene2'],'stock_gene1':s['gene1'],'diagnostic_gene1':r['gene1'],'only_gene1_order_differs':sorted(s['gene1'].split(','))==sorted(r['gene1'].split(','))})
out['characterized_differences']=differences
out['all_differences_only_gene1_list_order']=len(differences)==out['stock_only_row_instances']==out['diagnostic_only_row_instances'] and all(r['only_gene1_order_differs'] for r in differences)
out['interpretation']='Preserve distinct file hashes. If explicitly verified here, an intergenic alternative-gene-label ordering difference does not change coordinates, support, filters or other fields. This is not a byte-identical discarded file.'
(P/'diagnostics/discarded-output-comparison.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if not k.endswith('examples')})

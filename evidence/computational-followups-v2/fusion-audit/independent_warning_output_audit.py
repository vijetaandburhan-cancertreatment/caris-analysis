#!/usr/bin/env python3
"""Independent parsed-field multiset check; no source-call mutation."""
import collections,csv,hashlib,json,pathlib,resource,time
ws=pathlib.Path(__file__).resolve().parent;root=ws.parent;diag=root/'diagnostics';t=time.time()
paths=[root/'patient-D8/fusions.discarded.tsv',diag/'patient-D8-counter-rerun/fusions.discarded.tsv']
def digest(vals):return hashlib.blake2b(json.dumps(vals,ensure_ascii=True,separators=(',',':')).encode(),digest_size=32).digest()
def scan(path,collect=None):
 counts=collections.Counter();samples={};num=0
 with path.open() as f:
  reader=csv.reader(f,delimiter='\t');header=next(reader)
  for row in reader:
   assert len(row)==len(header);h=digest(row);num+=1
   if collect is None:counts[h]+=1
   elif h in collect:samples[h]=dict(zip(header,row))
 return header,num,counts,samples
h1,n1,c1,_=scan(paths[0]);h2,n2,c2,_=scan(paths[1]);assert h1==h2
missing1=c1-c2;missing2=c2-c1
_,_,_,samples1=scan(paths[0],missing1);_,_,_,samples2=scan(paths[1],missing2)
# Match changed rows on every field except gene1; diagnose, do not silently normalize new biology.
def k(row):return digest([(name,row[name]) for name in h1 if name!='#gene1'])
buckets=collections.defaultdict(list)
for h,row in samples2.items():buckets[k(row)].append((h,row))
differences=[];unmatched=[]
for h,left in samples1.items():
 cand=buckets[k(left)]
 if len(cand)!=1:unmatched.append({'side':'stock','row':left,'matching_diagnostic_rows_excluding_gene1':len(cand)});continue
 hh,right=cand.pop();changes={key:{'stock':left[key],'diagnostic':right[key]} for key in h1 if left[key]!=right[key]}
 names1=left['#gene1'].split(',');names2=right['#gene1'].split(',');equal_tokens=collections.Counter(names1)==collections.Counter(names2)
 differences.append({'changes':changes,'same_all_other_fields':set(changes)=={'#gene1'},'gene1_tokens_equal_including_multiplicity':equal_tokens,'stock_multiplicity':missing1[h],'diagnostic_multiplicity':missing2[hh],'breakpoint1':left['breakpoint1'],'breakpoint2':left['breakpoint2'],'gene2':left['gene2'],'filters':left['filters']})
for v in buckets.values():
 for _,row in v:unmatched.append({'side':'diagnostic','row':row})
passed=n1==n2 and not unmatched and all(d['same_all_other_fields'] and d['gene1_tokens_equal_including_multiplicity'] and d['stock_multiplicity']==d['diagnostic_multiplicity'] for d in differences)
w=json.load(open(diag/'warning-analysis.json'));reasons=w['reason_counts'];observed={};log=diag/'patient-D8-counter-rerun/Arriba.log'
for line in log.open():
 col=line.strip().split('\t')
 if col[0]=='ARRIBA_DIAG_REASON_TOTAL':observed[col[1]]=int(col[2])
assert all(observed.get(k,0)==v for k,v in reasons.items())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
accepted_equal=(root/'patient-D8/fusions.tsv').read_bytes()==(diag/'patient-D8-counter-rerun/fusions.tsv').read_bytes()
r={'status':'passed' if passed and accepted_equal else 'failed','method':'Separate CSV field parser with all-field BLAKE2b-256 row-multisets, retaining multiplicity and original read-name ordering. Every changed row independently matched on all fields except gene1; literal gene1 token multisets compared. No unconditional gene normalization.','discarded_rows':[n1,n2],'discarded_headers_equal':h1==h2,'all_field_identical_row_instances':n1-sum(missing1.values()),'stock_differing_row_instances':sum(missing1.values()),'diagnostic_differing_row_instances':sum(missing2.values()),'differences':differences,'unmatched':unmatched,'semantic_field_equality_except_equal_gene1_token_order':passed,'accepted_bytes_equal':accepted_equal,'accepted_sha256':sha(root/'patient-D8/fusions.tsv'),'stock_discarded_sha256':sha(paths[0]),'diagnostic_discarded_sha256':sha(paths[1]),'warning_counters_independently_parsed_from_log':observed,'reason_sum':sum(observed.values()),'warning_source_JSON_values_match':True,'expected_warning':w['stock_malformed_warning'],'units_review':'61,161 rejected paired alignment groups with group size other than2or3;5,303 groups where split/supplementary overlap could not be resolved;37,468 individual supplementary records clipped at the wrong end. Same fragment can contribute across records/hits/stages; total is not unique reads, molecules, loss percentage or sensitivity.','subsampling_review':'Counts report cap-branch encounters/distinct keys, not absolute missed-fragment counts. Filtering thresholds and default300cap unchanged.','limits':['Accepted byte identity and parsed-field discarded equivalence verify this counter rerun did not change the reported calls; they do not validate biological calls or eliminate sensitivity loss from filtering.','The absence of fatal BAM-decoder errors rules out that observed failure route, not every possible upstream sequencing/alignment problem.','Sparse D8 indexing remains separately non-equivalent to compact D1 control. No whole-input detection sensitivity was estimated.'],'elapsed_seconds':time.time()-t,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(ws/'independent-warning-output-audit.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));assert r['status']=='passed'

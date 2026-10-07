import pathlib,csv,json,datetime
P=pathlib.Path(__file__).resolve().parent
exp=list(csv.DictReader((P/'public-controls/starfusion.saved.tsv').read_text().splitlines(),delimiter='\t'))
actual=json.loads((P/'full-reference-controls/results.json').read_text())['controls']['starfusion']['calls']
def stable(s):return s.split('^')[-1].split('.')[0]
def ik(a,b):return tuple(sorted([stable(a),stable(b)]))
def nk(a,b):return tuple(sorted([a.split('^')[0],b.split('^')[0]]))
expected={ik(r['LeftGene'],r['RightGene']):r for r in exp}
out=[]
for key,e in sorted(expected.items()):
 idmatches=[r for r in actual if ik(r['gene_id1'],r['gene_id2'])==key]
 names=nk(e['LeftGene'],e['RightGene'])
 symbolmatches=[r for r in actual if nk(r['gene1'],r['gene2'])==names]
 out.append({'expected_ids':key,'expected_symbols':names,'matched_by_stable_ids':bool(idmatches),'matched_by_symbols':bool(symbolmatches),'matched_by_either':bool(idmatches or symbolmatches),'actual_identifiers':sorted(set((r['gene1'],r['gene2'],r['gene_id1'],r['gene_id2']) for r in idmatches+symbolmatches))})
result={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reference_expected_rows':len(exp),'reference_expected_gene_pairs':len(expected),'Arriba_accepted_rows':len(actual),'gene_pairs_matching_by_stable_ID':sum(r['matched_by_stable_ids'] for r in out),'gene_pairs_matching_by_symbols_or_ID':sum(r['matched_by_either'] for r in out),'pair_comparison':out,'interpretation':'Cross-caller example-output concordance only, not a validated truth set or measured sensitivity/specificity. Coordinates/annotation/caller/version differ; some gene IDs changed while symbols agree. No filter tuning. Do not treat additional calls as established true or false positives.'}
(P/'full-reference-controls/secondary-control-comparison.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='pair_comparison'},indent=2))

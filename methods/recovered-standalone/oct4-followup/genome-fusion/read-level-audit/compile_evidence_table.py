#!/usr/bin/env python3
import argparse,collections,csv,json,pathlib

def main():
 p=argparse.ArgumentParser();p.add_argument('--audit',required=True);p.add_argument('--reference',required=True);p.add_argument('--out',required=True);a=p.parse_args();audit=json.load(open(a.audit));ref=json.load(open(a.reference));dest=pathlib.Path(a.out);dest.mkdir(parents=True,exist_ok=True);lookup=collections.defaultdict(set)
 for i,patt in enumerate(ref['patterns']):
  for r in patt['roles']:lookup[(r['candidate'],r['kind'],r['arm'])].add(i)
 rows=[];details=[]
 for r in audit['rows']:
  cid=r['audit_key'];s=audit['statistics'][cid];o={'candidate_id':cid,'source_set':r['source_set'],'gene1':r['gene1'],'gene2':r['gene2'],'breakpoint1':r['breakpoint1'],'breakpoint2':r['breakpoint2'],'confidence':r['confidence'],'frame':r['reading_frame'],'type':r['type'],'caller_split_support':r['split_support'],'caller_discordant_mates':r['discordant_mates'],'caller_named_pairs':r['distinct_named_reads'],'caller_named_pairs_recovered':s['caller_named_pairs_recovered'],'caller_missing_names':len(s['caller_names_missing']),'exact_pattern_arms_available':','.join(str(x['arm_bases']) for x in r['junction_patterns']),'caller_filters':r['filters'],'caller_tags':r['tags']};rd={'candidate':r,'reference_matches':{}}
  for arm in (15,20,25):
   for q in (20,30):
    cs=s['raw_exact'][str(arm)][f'Q{q}']
    for field in ['pairs','caller_named_pairs','sequence_families','caller_named_sequence_families','orientation_normalized_sequence_families']:o[f'arm{arm}_Q{q}_{field}']=cs[field]
   for kind in ['junction','left_anchor','right_anchor']:
    ii=lookup[(cid,kind,arm)];gcount=sum(ref['genome']['counts'].get(str(i),0) for i in ii);txcount=sum(ref['transcripts']['matching_transcripts'].get(str(i),0) for i in ii);genes=sorted({g for i in ii for g in ref['transcripts']['matching_genes'].get(str(i),{})})
    o[f'{kind}_arm{arm}_genome_occurrences']=gcount;o[f'{kind}_arm{arm}_normal_transcript_pattern_occurrences']=txcount;o[f'{kind}_arm{arm}_normal_genes']=';'.join(genes)
    rd['reference_matches'][f'{kind}_arm{arm}']={'pattern_indices':sorted(ii),'genome_occurrences':gcount,'normal_transcript_pattern_occurrences':txcount,'normal_genes':genes,'genome_first_loci':[{'pattern_index':i,**h} for i in ii for h in ref['genome']['first_loci'].get(str(i),[])],'transcript_first_loci':[{'pattern_index':i,**h} for i in ii for h in ref['transcripts']['first_hits'].get(str(i),[])]}
  o['review_status']='unresolved_pending_scientific_review';o['status_reason']='No automatic biological/actionability inference from exact-match counts or gene names.';rows.append(o);details.append(rd)
 with (dest/'candidate-evidence.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['candidate_id'],delimiter='\t');w.writeheader();w.writerows(rows)
 (dest/'candidate-reference-details.json').write_text(json.dumps(details,indent=2));(dest/'candidate-evidence.json').write_text(json.dumps(rows,indent=2));print(json.dumps({'rows':len(rows),'status':'evidence_compiled_review_pending'}))
if __name__=='__main__':main()

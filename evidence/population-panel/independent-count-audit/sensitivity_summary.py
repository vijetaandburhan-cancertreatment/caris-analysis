from pathlib import Path
import csv,json,collections,hashlib,datetime
O=Path(__file__).resolve().parent
MODES=['baseline','clean_no_softclip_edge5','MAPQ60','NH_unique_if_tagged']
def load(p):return list(csv.DictReader(p.open(),delimiter='\t'))
d={a:{(r['chrom'],int(r['pos1'])):r for r in load(O/f'{a}-direct-counts.tsv')} for a in ['DNA','RNA']}
annotation={(r['chrom'],int(r['pos1'])):r for r in load(O.parent/'joint-callable.tsv')}
assert set(d['DNA'])==set(d['RNA'])==set(annotation)
summary={};events={};fixed={}
for mode in MODES:
 ctr=collections.Counter();ret=collections.defaultdict(collections.Counter);ev=[];fixedctr=collections.Counter()
 for k,dr in d['DNA'].items():
  rr=d['RNA'][k];x={a:{b:int(d[a][k][mode+'_'+b]) for b in 'ACGT'} for a in ['DNA','RNA']};dep={a:sum(x[a].values()) for a in x}
  ctr['baseline_joint_loci']+=1
  if min(dep.values())>=20:
   ctr['mode_joint_callable20']+=1;ref,alt=dr['ref'],dr['alt'];de=dep['DNA']
   cat='ref_dominant' if x['DNA'][ref]*100>=98*de else 'alt_dominant' if x['DNA'][alt]*100>=98*de else 'ref_alt_mixed' if min(x['DNA'][ref],x['DNA'][alt])*10>=de and (x['DNA'][ref]+x['DNA'][alt])*100>=95*de else 'residual';ret[cat]['n']+=1
   if cat.endswith('_dominant'):
    same=ref if cat.startswith('ref') else alt
    for cut in [90,98]:ret[cat]['retain'+str(cut)]+=x['RNA'][same]*100>=cut*dep['RNA']
   oldcat=annotation[k]['DNA_category']
   if oldcat.endswith('_dominant'):
    fixedctr['baseline_dominant_still_joint20']+=1;same=ref if oldcat.startswith('ref') else alt
    for cut in [90,98]:fixedctr['RNA_retains_original_DNA_dominant_at'+str(cut)]+=x['RNA'][same]*100>=cut*dep['RNA']
  if min(dep.values())>=50:
   for b in 'ACGT':
    if x['DNA'][b]*100<=2*dep['DNA'] and x['RNA'][b]>=10 and x['RNA'][b]*10>=dep['RNA']:
     ev.append({'chrom':k[0],'pos1':k[1],'gene_symbols':annotation[k]['overlap_gene_symbols'],'unexpected_base':b,'DNA':x['DNA'],'RNA':x['RNA'],'DNA_depth':dep['DNA'],'RNA_depth':dep['RNA'],'RNA_fraction':x['RNA'][b]/dep['RNA'],'DNA_fraction':x['DNA'][b]/dep['DNA']})
 summary[mode]={'attrition':dict(ctr),'categories':{k:dict(v) for k,v in ret.items()},'exact_strong_allele_events':len(ev),'exact_strong_loci':len(set((e['chrom'],e['pos1']) for e in ev))};events[mode]=ev;fixed[mode]=dict(fixedctr)
# Preserve all baseline events with every sensitivity, even when thresholds no longer met.
details=[]
for e in events['baseline']:
 k=(e['chrom'],e['pos1']);z=dict(e);z['sensitivity_all_counts']={}
 for mode in MODES[1:]:
  z['sensitivity_all_counts'][mode]={a:{b:int(d[a][k][mode+'_'+b]) for b in 'ACGT'} for a in ['DNA','RNA']}
 details.append(z)
result={'status':'PASS','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'method':'Exact per-base allele screen over all4015 baseline jointly callable sites, including mixed/residual DNA. Both ACGT depths>=50; each unexpected base DNA<=2%, RNA>=10% and>=10QNAME counts. Primary aggregation is not used for exact-base triage. Mode summaries recategorize only loci retaining depth20; fixed-baseline summaries retain original dominant-allele designation.','mode_summary':summary,'fixed_baseline_dominant_sensitivity':fixed,'exact_events_by_mode':events,'baseline_events_all_sensitivities':details,'input_hashes':{a:hashlib.file_digest((O/f'{a}-direct-counts.tsv').open('rb'),'sha256').hexdigest() for a in ['DNA','RNA']},'limitations':['Sensitivity filters do not establish gene of origin or explain discrepancies.','Drop below depth50 is loss of screen assessability, not resolution. For example APOD and EXOC3 remain discordant after clean filters but fail the depth50 screen.','Separate neighboring PABPC1 SNPs may arise from shared RNA fragments and are not independent evidence.','MAPQ60 and NH-unique results can remain unchanged for STAR alignments scored255; these tags do not establish unique molecular origin.','No identity probability, contamination fraction, clinical diagnosis or treatment inference.']}
(O/'sensitivity-and-exact-discordances.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'mode_summary':summary,'fixed_baseline':fixed},indent=2))

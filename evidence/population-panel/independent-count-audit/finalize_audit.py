from pathlib import Path
import csv,json,collections,hashlib,datetime
O=Path(__file__).resolve().parent;R=O.parent
MODES=['baseline','clean_no_softclip_edge5','MAPQ60','NH_unique_if_tagged']
def load(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def key(r):return r['chrom'],int(r['pos1'])
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
primary={key(r):r for r in load(R/'joint-callable.tsv')}
direct={a:{key(r):r for r in load(O/f'{a}-direct-counts.tsv')} for a in ['DNA','RNA']}
def summarize(k,classification,unexpected_base=''):
 p=primary[k];out={f:p[f] for f in ['chrom','pos1','id','ref','alt','overlap_gene_symbols','gene_overlap_category','unique_gene_id','DNA_category']}
 out.update({'review_class':classification,'unexpected_base':unexpected_base,'baseline_RNA_same_DNA_dominant_fraction':p['RNA_same_allele_fraction'],'both_baseline_depths_ge50':int(min(int(p['DNA_ACGT_depth']),int(p['RNA_ACGT_depth']))>=50)})
 for mode in MODES:
  for a in ['DNA','RNA']:
   r=direct[a][k]
   out[f'{mode}_{a}_ACGT']=','.join(str(r[mode+'_'+b]) for b in 'ACGT')
   out[f'{mode}_{a}_ACGT_depth']=r[mode+'_ACGT_depth']
   out[f'{mode}_{a}_conflicting_names']=r[mode+'_fragment_discordant']
 return out
ex=[]
for k,p in primary.items():
 if p['DNA_category'].endswith('_dominant') and p['retain90']=='0':
  strong=p['triage_candidate']=='1';ex.append(summarize(k,'strong_baseline_screen' if strong else 'lower_RNA_depth_exception'))
assert len(ex)==13 and sum(x['review_class']=='lower_RNA_depth_exception' for x in ex)==7
write(O/'all13-below90-dominant-sites.tsv',ex)
ev=json.loads((O/'sensitivity-and-exact-discordances.json').read_text());baseline={key(r) for r in ev['exact_events_by_mode']['baseline']};mq=[]
for e in ev['exact_events_by_mode']['MAPQ60']:
 k=key(e)
 if k not in baseline:mq.append(summarize(k,'MAPQ60_only_not_baseline_flag',e['unexpected_base']))
assert len(mq)==3;write(O/'three-MAPQ60-only-flags.tsv',mq)
# Independently intersect the two retained panel-chromosome genes using separately
# GENCODE-derived gene feature inventory; no target BAM or primary interval code used.
inv=load(R/'public-annotation/overlapping-gene-inventory.tsv')
regions={r['gene_symbol']:r for r in inv if r['gene_symbol'] in ['BAP1','RASA1']}
allrows=load(R/'all-public-loci-counts.tsv');retained=[]
for r in allrows:
 for gene,g in regions.items():
  if r['chrom']==g['chrom'] and int(g['start1'])<=int(r['pos1'])<=int(g['end1']):retained.append((gene,r))
intersection=[r for gene,r in retained if int(r['joint_callable'])]
check={'status':'PASS','method':'Direct inclusive containment using GENCODE-derived BAP1/RASA1 gene spans and completed full-panel count table, without importing the other agent interval function or opening any BAM.','public_loci':len(retained),'by_gene':dict(collections.Counter(g for g,r in retained)),'joint_loci':len(intersection),'DNA_callable20':sum(int(r['DNA_ACGT_depth'])>=20 for g,r in retained),'RNA_callable20':sum(int(r['RNA_ACGT_depth'])>=20 for g,r in retained),'max_RNA_ACGT_depth':max(int(r['RNA_ACGT_depth']) for g,r in retained),'regions':regions,'interpretation':'No eligible cross-alignment comparison; zero rows cannot be interpreted as concordance.'}
assert (check['public_loci'],check['joint_loci'],check['DNA_callable20'],check['RNA_callable20'],check['max_RNA_ACGT_depth'])==(73,0,3,0,9)
(O/'independent-target-intersection.json').write_text(json.dumps(check,indent=2)+'\n')
schema={'count_string_order':'A,C,G,T, comma-separated counts of unambiguous qualifying QNAME units; not UMI molecules.','baseline_filters':'MAPQ>=30,BQ>=25,proper pair,exclude0xF0C; conflicting ACGT mate-name observations discarded.','scope':'All13 DNA-dominant primary sites below90% RNA same-allele retention, plus separate three MAPQ60-only flags. Baseline sites are not replaced by sensitivity results.','lower_depth':'All seven lower-depth exceptions have RNA ACGT depth24,31,26,39,26,22 or34 (under50). They are observable discrepancies but below the predeclared strong-screen depth floor, not explained artifacts.','observed_row_vs_ACGT':'A primary pileup row can contain deletion/refskip/nonACGT evidence or insufficient qualifying ACGT counts. Observed-row totals (DNA62358,RNA78071) are not callable-site counts (DNA35885,RNA4072 at ACGT depth20); only4015 qualify in both.','MAPQ_sensitivity':'MAPQ60 can remove one DNA allele while leaving STAR RNA unchanged, producing new fraction flags. This is a mapping/selection sensitivity, not evidence of a newly discovered biological variant.','new_alignment':'No jointly callable public-panel sites intersect retained BAP1/RASA1 target intervals, so independent alignment replication is not evaluable.'}
(O/'exception-tables-schema.json').write_text(json.dumps(schema,indent=2)+'\n')
p=O/'findings.txt';text=p.read_text().split('\nComplete exception accounting:')[0];text+='\nComplete exception accounting: all13-below90-dominant-sites.tsv contains all13 primary dominant sites below90% with exact A/C/G/T, depth, conflicts and every sensitivity. The seven lower-RNA-depth exceptions are WWC1 (24), RGS14 (31), PTCH1 (26), NOTCH1 (39), MYO18A (26), NF1 (22), and SEPTIN9 (34); parenthesized values are baseline RNA ACGT depth. They are retained as observations, not dismissed. three-MAPQ60-only-flags.tsv separately preserves all three sensitivity-only sites. Count-string order is A,C,G,T, documented in exception-tables-schema.json.\n\nObserved rows do not equal callable ACGT coverage: DNA has62,358 observed pileup rows and35,885 depth20 sites; RNA has78,071 observed rows and4,072 depth20 sites. Rows with low/non-ACGT evidence remain in attrition. An independent interval check confirms73 public-panel sites intersect the retained BAP1/RASA1 alignment intervals, but none is RNA-depth20 or jointly callable (maximum RNA depth9). Consequently the second-alignment check is not evaluable, not concordant.\n';p.write_text(text)
coverage={}
for assay in ['DNA','RNA']:
    rs=[r for ch in ['chr2','chr3','chr5','chr8','chr9','chr17'] for r in load(R/assay/f'{ch}.tsv')]
    depth=[sum(int(r['fragment_'+v]) for v in ['ref','alt','other']) for r in rs]
    coverage[assay]={'observed_pileup_rows':len(rs),'at_least1_ACGT_QNAME':sum(x>=1 for x in depth),'at_least20_ACGT_QNAME':sum(x>=20 for x in depth),'zero_ACGT_rows':sum(x==0 for x in depth)}
assert coverage['DNA']['at_least1_ACGT_QNAME']==62262 and coverage['RNA']['at_least1_ACGT_QNAME']==25707
(O/'independent-observed-versus-ACGT-coverage.json').write_text(json.dumps({'status':'PASS','method':'Direct summation of fragment_ref+fragment_alt+fragment_other in each original per-assay chromosome TSV, without treating observed rows or nonACGT/refskip symbols as ACGT coverage.','counts':coverage},indent=2)+'\n')
with p.open('a') as f:f.write('\nExact nonzero-base coverage from the per-assay TSVs: DNA62,262 loci have>=1ACGT QNAME versus62,358 observed rows; RNA25,707 have>=1ACGT QNAME versus78,071 observed rows. Thus96 DNA and52,364 RNA observed rows have zero unambiguous ACGT QNAME counts. The independent 4,015-site count comparison intentionally excludes symbolic skip/deletion non-ACGT accounting, consistent with the ACGT denominator.\n')
files=[x for x in O.iterdir() if x.is_file() and x.name!='independent-audit-manifest.json']
man={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'design_sha256':'90a153f630375103363792037618d93a864b01bad029575f36d6a3a382ee0bfa','status':'PASS exact count and summary QC; origin interpretation separate','outputs':{x.name:{'bytes':x.stat().st_size,'sha256':hashlib.file_digest(x.open('rb'),'sha256').hexdigest()} for x in sorted(files)}}
(O/'independent-audit-manifest.json').write_text(json.dumps(man,indent=2)+'\n')
print({'all13':len(ex),'MAPQ60_only':len(mq),'target_intersection':check['joint_loci']})

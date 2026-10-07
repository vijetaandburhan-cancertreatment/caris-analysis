"""Independent arithmetic on frozen public universe and primary per-assay counts."""
from pathlib import Path
import collections,csv,json,hashlib,statistics,datetime
O=Path(__file__).resolve().parent;R=O.parent
CH=['chr2','chr3','chr5','chr8','chr9','chr17']
def load(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def key(r):return (r['chrom'],int(r['pos1']))
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
panel=load(R/'input/public-common-SNP-panel.tsv');panel.sort(key=lambda r:(CH.index(r['chrom']),int(r['pos1'])))
counts={a:{key(r):r for c in CH for r in load(R/a/f'{c}.tsv')} for a in ['DNA','RNA']}
annot={key(r):r for r in load(R/'public-annotation/public-panel-gene-annotations.tsv')}
primary={key(r):r for r in load(R/'all-public-loci-counts.tsv')}
assert len(panel)==len(primary)==165782
cats=collections.Counter();catstats=collections.defaultdict(collections.Counter);spstats=collections.defaultdict(collections.Counter)
last={};spkeys=set();joint=set();dominant=set();dis=[];attr=collections.Counter();chrattr=collections.defaultdict(collections.Counter)
for p in panel:
 k=key(p);ch,po=k;sr=primary[k];attr['selected']+=1;chrattr[ch]['selected']+=1
 if ch not in last or po-last[ch]>=50000:spkeys.add(k);last[ch]=po
 assert bool(int(sr['frozen_spaced50kb']))==(k in spkeys)
 x={a:{b:int(counts[a].get(k,{}).get('fragment_'+b,0)) for b in ['ref','alt','other']} for a in ['DNA','RNA']}
 d={a:sum(x[a].values()) for a in x}
 for a in x:
  attr[a+'_observed']+=k in counts[a];attr[a+'_callable']+=d[a]>=20;chrattr[ch][a+'_observed']+=k in counts[a];chrattr[ch][a+'_callable']+=d[a]>=20
 cat='not_joint_callable'
 if min(d.values())>=20:
  joint.add(k);attr['joint']+=1;chrattr[ch]['joint']+=1
  ref,alt=x['DNA']['ref'],x['DNA']['alt'];dep=d['DNA']
  cat='ref_dominant' if ref*100>=98*dep else 'alt_dominant' if alt*100>=98*dep else 'ref_alt_mixed' if min(ref,alt)*10>=dep and (ref+alt)*100>=95*dep else 'residual'
  cats[cat]+=1;catstats[cat]['n']+=1
  if k in spkeys:spstats[cat]['n']+=1
  if cat.endswith('_dominant'):
   dominant.add(k);same='ref' if cat.startswith('ref') else 'alt';opposite='alt' if same=='ref' else 'ref'
   for cut in [90,98]:
    val=x['RNA'][same]*100>=cut*d['RNA'];catstats[cat]['retain'+str(cut)]+=val
    if k in spkeys:spstats[cat]['retain'+str(cut)]+=val
    if val!=bool(int(sr['retain'+str(cut)])):dis.append({'site':k,'field':'retain'+str(cut),'integer':val,'primary':sr['retain'+str(cut)]})
   val=x['RNA'][opposite]*100>=98*d['RNA'];catstats[cat]['opposite98']+=val
   if k in spkeys:spstats[cat]['opposite98']+=val
 if cat!=sr['DNA_category']:dis.append({'site':k,'field':'DNA_category','integer':cat,'primary':sr['DNA_category']})
def breadth(keys):
 by=collections.Counter(c for c,p in keys);genes={annot[k]['unique_gene_id'] for k in keys if annot[k]['unique_gene_id']!='.'};bins={(c,(p-1)//1000000) for c,p in keys};gaps=[]
 for c in CH:
  ps=sorted(p for ch,p in keys if ch==c);gaps.extend(b-a for a,b in zip(ps,ps[1:]))
 return {'loci':len(keys),'chromosomes':dict(by),'unique_genes':len(genes),'oneMb_bins':len(bins),'distance_n':len(gaps),'distance_min':min(gaps) if gaps else None,'distance_median':statistics.median(gaps) if gaps else None,'distance_max':max(gaps) if gaps else None,'distance_under1000':sum(x<1000 for x in gaps)}
old={key(r):r for r in load(R/'input/public-common-SNP-alleles.annotated.tsv')}
cols=['read_ref','read_alt','fragment_ref','fragment_alt','fragment_other','fragment_discordant','fragment_nonACGT','fragment_depth_refalt','alt_forward_reads','alt_reverse_reads','ref_forward_reads','ref_reverse_reads','fragment_ref_BP5_gt5','fragment_alt_BP5_gt5','reported_pileup_depth']
diff=[]
for k,r in old.items():
 if k not in counts['DNA']:diff.append({'site':k,'missing':True});continue
 for field in cols:
  if int(r[field])!=int(counts['DNA'][k][field]):diff.append({'site':k,'field':field,'old':r[field],'new':counts['DNA'][k][field]})
result={'status':'PASS' if not dis else 'REVIEW_REQUIRED','method':'Reconstruct from frozen public panel and per-assay TSVs; exact integer threshold arithmetic, including ref+alt>=95%, independent of primary merged categorical labels. Compare all merged categories, specified allele-retention flags and preselected spacing.','panel_loci':len(panel),'joint_loci':len(joint),'categories':dict(cats),'category_retention':{k:dict(v) for k,v in catstats.items()},'attrition':dict(attr),'by_chr_attrition':{k:dict(v) for k,v in chrattr.items()},'joint_breadth':breadth(joint),'dominant_breadth':breadth(dominant),'spaced_public':len(spkeys),'spaced_joint_breadth':breadth(spkeys&joint),'spaced_retention':{k:dict(v) for k,v in spstats.items()},'summary_differences':dis,'old_DNA_comparison':{'old_loci':len(old),'new_loci':len(counts['DNA']),'new_sites':len(set(counts['DNA'])-set(old)),'old_missing':len(set(old)-set(counts['DNA'])),'integer_fields_per_locus':len(cols),'differences':diff},'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
save(O/'independent-summary-audit.json',result)
print(json.dumps({k:v for k,v in result.items() if k not in ['by_chr_attrition','old_DNA_comparison']},indent=2))

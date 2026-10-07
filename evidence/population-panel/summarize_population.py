import pathlib,csv,json,collections,statistics,hashlib
O=pathlib.Path(__file__).resolve().parent
W=pathlib.Path('/Users/burhanazeem/Documents/Codex/2026-09-05/finances-plugin-finances-openai-curated-remote-3')
C=['chr2','chr3','chr5','chr8','chr9','chr17']
def load(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def key(r):return r['chrom'],int(r['pos1'])
def write(p,rows):
 if not rows:p.write_text('');return
 with p.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def js(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert json.loads((O/'count-complete.json').read_text())['status']=='PASS'
panel=load(O/'public-annotation/public-panel-gene-annotations.tsv');panel.sort(key=lambda r:(C.index(r['chrom']),int(r['pos1'])))
counts={a:{key(r):r for ch in C for r in load(O/a/f'{ch}.tsv')} for a in ['DNA','RNA']}
old={key(r):r for r in load(O/'input/public-common-SNP-alleles.annotated.tsv')}
cols=['read_ref','read_alt','fragment_ref','fragment_alt','fragment_other','fragment_discordant','fragment_nonACGT','fragment_depth_refalt','alt_forward_reads','alt_reverse_reads','ref_forward_reads','ref_reverse_reads','fragment_ref_BP5_gt5','fragment_alt_BP5_gt5','reported_pileup_depth']
diff=[]
for k,r in old.items():
 n=counts['DNA'].get(k)
 if n is None:diff.append({'chrom':k[0],'pos1':k[1],'field':'missing','old':'present','new':'absent'});continue
 for c in cols:
  if int(r[c])!=int(n[c]):diff.append({'chrom':k[0],'pos1':k[1],'field':c,'old':r[c],'new':n[c]})
newkeys=set(counts['DNA'])-set(old)
write(O/'old-DNA-count-differences.tsv',diff)
js(O/'old-DNA-comparison.json',{'prior_observed_loci':len(old),'new_observed_loci':len(counts['DNA']),'integer_fields_per_row':len(cols),'compared_integer_fields':len(old)*len(cols),'differences':len(diff),'new_loci':len(newkeys),'old_missing_loci':len(set(old)-set(counts['DNA'])),'old_maximum_reported_postquality_depth':max(int(x['reported_pileup_depth']) for x in old.values()),'new_maximum_reported_postquality_depth':max(int(x['reported_pileup_depth']) for x in counts['DNA'].values()),'interpretation':'Exact count equality, if observed, supports no depth-cap effect for these particular final fields; this is an empirical rerun check, not inference from postquality depth alone.'})
old335={key(r) for r in load(W/'work/oct1-deep/dna-rna-consistency/allele-consistency.tsv')}
last={};joined=[];triage=[]
for r0 in panel:
 r=dict(r0);k=key(r);ch,pos=k
 r['frozen_spaced50kb']=int(pos-last.get(ch,-50000)>=50000)
 if r['frozen_spaced50kb']:last[ch]=pos
 r['prior335_overlap']=int(k in old335)
 for a in ['DNA','RNA']:
  n=counts[a].get(k,{})
  r[a+'_observed_row']=int(bool(n))
  for c in cols:r[a+'_'+c]=int(n.get(c,0))
  d=sum(r[a+'_fragment_'+c] for c in ['ref','alt','other']);r[a+'_ACGT_depth']=d
  for c in ['ref','alt','other']:r[a+'_'+c+'_fraction']=r[a+'_fragment_'+c]/d if d else ''
  r[a+'_callable20']=int(d>=20)
 r['joint_callable']=int(r['DNA_callable20'] and r['RNA_callable20']);r['DNA_category']='not_joint_callable';r['RNA_same_allele_fraction']='';r['retain90']='';r['retain98']='';r['opposite98']=0;r['triage_candidate']=0;r['third_allele_excess']=0
 if r['joint_callable']:
  dr,da,do=[r['DNA_'+c+'_fraction'] for c in ['ref','alt','other']];rr,ra,ro=[r['RNA_'+c+'_fraction'] for c in ['ref','alt','other']]
  cat='ref_dominant' if dr>=.98 else 'alt_dominant' if da>=.98 else 'ref_alt_mixed' if min(dr,da)>=.10 and dr+da>=.95 else 'residual';r['DNA_category']=cat
  if cat in ['ref_dominant','alt_dominant']:
   c='ref' if cat=='ref_dominant' else 'alt';other='alt' if c=='ref' else 'ref';same=r['RNA_'+c+'_fraction'];unexpected=r['RNA_ACGT_depth']-r['RNA_fragment_'+c]
   r['RNA_same_allele_fraction']=same;r['retain90']=int(same>=.90);r['retain98']=int(same>=.98);r['opposite98']=int(r['RNA_'+other+'_fraction']>=.98)
   r['triage_candidate']=int(min(r['DNA_ACGT_depth'],r['RNA_ACGT_depth'])>=50 and unexpected>=10 and unexpected/r['RNA_ACGT_depth']>=.10)
  r['third_allele_excess']=int(min(r['DNA_ACGT_depth'],r['RNA_ACGT_depth'])>=50 and do<=.02 and ro>=.10 and r['RNA_fragment_other']>=10)
  if r['triage_candidate'] or r['third_allele_excess']:triage.append(r)
 joined.append(r)
joint=[r for r in joined if r['joint_callable']];dom=[r for r in joint if r['DNA_category'].endswith('_dominant')]
def distribution(rs):
 uniq=collections.Counter(r['unique_gene_id'] for r in rs if r['unique_gene_id']!='.')
 gaps=[]
 for ch in C:
  p=sorted(int(r['pos1']) for r in rs if r['chrom']==ch);gaps.extend(b-a for a,b in zip(p,p[1:]))
 return {'loci':len(rs),'chromosomes':dict(collections.Counter(r['chrom'] for r in rs)),'unique_1Mb_bins':len({(r['chrom'],(int(r['pos1'])-1)//1000000) for r in rs}),'uniquely_assigned_gene_count':len(uniq),'gene_overlap_categories':dict(collections.Counter(r['gene_overlap_category'] for r in rs)),'top_unique_genes':uniq.most_common(10),'scope':dict(collections.Counter(r['scope'] for r in rs)),'same_chromosome_consecutive_distance_bp':{'n':len(gaps),'min':min(gaps) if gaps else None,'median':statistics.median(gaps) if gaps else None,'max':max(gaps) if gaps else None,'under1000':sum(x<1000 for x in gaps)},'prior335_overlap':sum(r['prior335_overlap'] for r in rs)}
def retention(rs):
 out={}
 for cat in ['ref_dominant','alt_dominant','ref_alt_mixed','residual']:
  a=[r for r in rs if r['DNA_category']==cat];v={'n':len(a)}
  if cat.endswith('_dominant'):v.update({'retain90':sum(r['retain90'] for r in a),'retain98':sum(r['retain98'] for r in a),'opposite98':sum(r['opposite98'] for r in a),'triage_candidates':sum(r['triage_candidate'] for r in a)})
  else:v.update({'RNA_ref98':sum(r['RNA_ref_fraction']>=.98 for r in a),'RNA_alt98':sum(r['RNA_alt_fraction']>=.98 for r in a),'RNA_mixed10':sum(min(r['RNA_ref_fraction'],r['RNA_alt_fraction'])>=.10 for r in a)})
  out[cat]=v
 return out
attr=[]
for label,subset in [('ALL',joined)]+[(ch,[r for r in joined if r['chrom']==ch]) for ch in C]+[('scope:'+scope,[r for r in joined if r['scope']==scope]) for scope in sorted({r['scope'] for r in joined})]:
 a={'group':label,'selected':len(subset)}
 for assay in ['DNA','RNA']:
  for c in ['observed_row','callable20']:a[assay+'_'+c]=sum(r[assay+'_'+c] for r in subset)
 a['joint_callable']=sum(r['joint_callable'] for r in subset);a['joint_dominant']=sum(r['DNA_category'].endswith('_dominant') for r in subset);attr.append(a)
gd=distribution(dom);rr=retention(joint)
gate={'dominant_loci':len(dom)>=100,'ref_dominant':rr['ref_dominant']['n']>=20,'alt_dominant':rr['alt_dominant']['n']>=20,'chromosomes':len(gd['chromosomes'])>=4,'bins1Mb':gd['unique_1Mb_bins']>=20,'unique_genes':gd['uniquely_assigned_gene_count']>=30}
spaced=[r for r in joined if r['frozen_spaced50kb']];sj=[r for r in spaced if r['joint_callable']]
summary={'status':'PRIMARY_COUNTS_COMPLETE_PENDING_INDEPENDENT_AUDIT','design_sha256':sha(O/'design-frozen.json'),'full_panel_distribution':distribution(joined),'joint_distribution':distribution(joint),'dominant_distribution':gd,'attrition':attr,'joint_categories':rr,'operational_gate':gate,'gate_pass':all(gate.values()),'frozen_spaced50kb':{'all_selected':len(spaced),'joint_distribution':distribution(sj),'categories':retention(sj)},'triage_candidate_loci':len(triage),'opposite98_loci':sum(r['opposite98'] for r in joint),'third_allele_excess_loci':sum(r['third_allele_excess'] for r in joint),'triage_note':'Dominant-DNA candidates use aggregate unexpected RNA alleles as a sensitive screen; independent per-base audit must determine exact unexpected alleles. DNA-mixed RNA shifts do not establish incompatibility. More than20 candidates triggers fixed20 audit cap and limited unresolved reporting.'}
write(O/'all-public-loci-counts.tsv',joined);write(O/'joint-callable.tsv',joint);write(O/'discordance-triage.tsv',triage);write(O/'coverage-attrition.tsv',attr);js(O/'primary-summary.json',summary)
print(json.dumps({k:summary[k] for k in ['joint_categories','operational_gate','triage_candidate_loci','opposite98_loci','third_allele_excess_loci','frozen_spaced50kb']},indent=2))

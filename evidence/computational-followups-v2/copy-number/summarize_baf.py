"""Tumor-only common-SNP allelic imbalance: descriptive, not LOH calling."""
import pathlib,csv,json,gzip,collections,math
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=pathlib.Path(__file__).resolve().parent
def tsv(p,rows):
 if rows:
  with p.open('w') as f:w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
cen={}
with gzip.open(OUT/'public/cytoBandIdeo.txt.gz','rt') as f:
 for ln in f:
  c,s,e,b,st=ln.strip().split('\t')
  if st=='acen':cen.setdefault(c,[]).extend([int(s),int(e)])
cen={c:(min(z),max(z)) for c,z in cen.items()}
def arm(c,p):
 a,b=cen[c]
 return c+'p' if p<a else c+'q' if p>=b else 'centromere'
rows=[]
for p in sorted(OUT.glob('common-SNP-alleles.chr*.tsv')):
 for r in csv.DictReader(p.open(),delimiter='\t'):
  for k in ['pos1','read_ref','read_alt','fragment_ref','fragment_alt','fragment_other','fragment_discordant','fragment_nonACGT','fragment_depth_refalt','alt_forward_reads','alt_reverse_reads','ref_forward_reads','ref_reverse_reads','fragment_ref_BP5_gt5','fragment_alt_BP5_gt5','reported_pileup_depth']:r[k]=int(r[k])
  for k in ['fragment_ALT_fraction','minor_allele_fraction','public_AF','public_SAS_AF']:r[k]=float(r[k]) if r[k] else None
  r['arm']=arm(r['chrom'],r['pos1']-1);minor=min(r['fragment_alt'],r['fragment_ref']);other=r['fragment_other'];n=r['fragment_depth_refalt'];maf=r['minor_allele_fraction'];discord=r['fragment_discordant'];orientation=min(r['alt_forward_reads'],r['alt_reverse_reads'],r['ref_forward_reads'],r['ref_reverse_reads'])
  base=n>=30 and minor>=5 and maf>=.05 and other/max(n+other,1)<=.02 and discord/max(n+discord,1)<=.05 and orientation>=1
  strict=base and n>=50 and minor>=10 and orientation>=2 and min(r['fragment_ref_BP5_gt5']/max(r['fragment_ref'],1),r['fragment_alt_BP5_gt5']/max(r['fragment_alt'],1))>=.9
  r['candidate_heterozygous_base']=int(base);r['candidate_heterozygous_stringent']=int(strict);rows.append(r)
tsv(OUT/'public-common-SNP-alleles.annotated.tsv',rows)
rng=np.random.default_rng(20261004)
def stats(rs,label,flt):
 selected=[r for r in rs if r[flt]];v=np.array([r['minor_allele_fraction'] for r in selected]);blocks=collections.defaultdict(list)
 for r in selected:blocks[(r['chrom'],(r['pos1']-1)//1000000)].append(r['minor_allele_fraction'])
 bs=[float(np.median(v)) for v in blocks.values()]
 boot=[]
 if len(bs)>=5:
  a=np.array(bs)
  for k in range(3000):boot.append(float(np.median(rng.choice(a,len(a),replace=True))))
 return {'region':label,'filter':flt,'public_panel_sites_with_any_pileup_output':len(rs),'candidate_heterozygous_sites':len(v),'oneMb_blocks':len(bs),'site_MAF_median':float(np.median(v)) if len(v) else None,'site_MAF_Q25':float(np.quantile(v,.25)) if len(v) else None,'site_MAF_Q75':float(np.quantile(v,.75)) if len(v) else None,'block_median_MAF':float(np.median(bs)) if bs else None,'block_bootstrap95_low':float(np.quantile(boot,.025)) if boot else None,'block_bootstrap95_high':float(np.quantile(boot,.975)) if boot else None,'median_fragment_depth':float(np.median([r['fragment_depth_refalt'] for r in selected])) if selected else None}
summ=[]
for a in sorted(set(r['arm'] for r in rows)):
 for flt in ['candidate_heterozygous_base','candidate_heterozygous_stringent']:summ.append(stats([r for r in rows if r['arm']==a],a,flt))
ann=json.loads((OUT/'annotation.json').read_text())
for gene in ['BAP1','RASA1','MTAP','CDKN2A','CDKN2B','APC']:
 g=ann['genes'][gene]
 for pad in [0,250000,1000000,5000000]:
  rs=[r for r in rows if r['chrom']==g['chrom'] and g['start0']-pad<=r['pos1']-1<g['end0']+pad]
  for flt in ['candidate_heterozygous_base','candidate_heterozygous_stringent']:summ.append(stats(rs,f'{gene}_plusminus{pad}bp',flt))
tsv(OUT/'regional-allelic-imbalance-summary.tsv',summ)
examples=[]
for row in summ:
 if row['filter']!='candidate_heterozygous_stringent' or row['region'] not in ['chr3p','chr9p','chr9q','BAP1_plusminus1000000bp','RASA1_plusminus1000000bp']:continue
 m=row['site_MAF_median']
 if m is None:continue
 for total,minor in [(1,0),(2,0),(3,0),(3,1),(4,0),(4,1),(5,1)]:
  # m=(1-p+p*minor)/(2-2p+p*total).
  den=m*(total-2)-(minor-1);p=(1-2*m)/den if den else None
  if p is not None and 0<=p<=1:
   bulk=2*(1-p)+p*total
   examples.append({'region':row['region'],'observed_median_minor_fraction':m,'assumed_clonal_tumor_total_copies':total,'assumed_tumor_minor_copies':minor,'implied_tumor_fraction_under_model':p,'implied_bulk_total_copies_under_model':bulk,'expected_depth_relative_to_true_diploid_normal':bulk/2,'interpretation':'illustrative mathematically compatible state; not an inferred copy call or purity estimate'})
tsv(OUT/'allelic-imbalance-nonidentifiability-examples.tsv',examples)
(OUT/'allelic-imbalance-method.json').write_text(json.dumps({'public_selection':'independent population-common panel, globalAF0.05–0.95','base_filter':'ref+altfragments>=30;bothalleles>=5;minorAF>=0.05;otherAF<=2%;discordant<=5%;bothorientationsforbothalleles>=1','stringent_filter':'basefilter plusdepth>=50;bothalleles>=10;bothorientations>=2;>=90%fragmentsfor eachallelehaveatleastonebase>5ntfromread5primeend','sampling':'No matchednormal. Sites with bothalleles become candidateheterozygotes, not confirmed constitutional heterozygotes. Complete LOH in puretumor could disappear from this selectedsubset.','CI':'3000 fixedseed bootstrap of1Mb blockmedians; descriptive samplingstability only, not confidenceinLOH or adjusted for technicalcapture bias. Regions<5blocks have no CI.','purityexamples':'Assume germlineAB normaldiploid, one homogeneous tumorstate, no mappingbias; m=[(1-p)+p*b]/[2(1-p)+p*C]. Multiple(C,b,p) produce same m. Bulkdepthratios shown against a TRUE normal, which is absent.','caveats':['PopulationLD and correlatedcapture sites prevent treating everymarker as an independentreplicate.','FFPE/mapping/reference/capture allele biases and subclonality are not eliminated.','Markers not physically phased to driver variants.','Arm-median estimates may average distinct segments.','No tumor-only allele counts prove a BAP1 secondhit, somaticstatus,orMTAPhomozygousloss.']},indent=2)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'font.size':10})
fig,axes=plt.subplots(4,1,figsize=(12,11),sharey=True)
for ax,ch in zip(axes,['chr3','chr5','chr9','chr17']):
 rs=[r for r in rows if r['chrom']==ch and r['candidate_heterozygous_stringent']];x=[r['pos1']/1e6 for r in rs];y=[r['fragment_ALT_fraction'] for r in rs];ax.scatter(x,y,s=8,color='#3c6975',alpha=.5)
 ax.axhline(.5,color='#444444',lw=.8,ls=':');a,b=cen[ch];ax.axvspan(a/1e6,b/1e6,color='#dddddd',alpha=.6);ax.set_ylim(-.02,1.02);ax.set_title(f'{ch} — {len(rs):,} candidate heterozygous common SNPs',loc='left',weight='bold');ax.set_ylabel('Alternate-allele fraction');ax.set_xlabel('Genomic position (Mb, hg38)')
 for gene,pos in [('BAP1',52.40855),('RASA1',87.33256),('APC',112.843324),('MTAP / CDKN2A/B',21.9)]:
  if (ch=='chr3' and gene=='BAP1') or (ch=='chr5' and gene in ['RASA1','APC']) or (ch=='chr9' and gene.startswith('MTAP')):ax.axvline(pos,color='#111111',ls='--',lw=.9);ax.text(pos+.8,.88,gene,fontsize=9)
fig.suptitle('Independent common-SNP allele balance in tumor DNA\nQuery-name fragment fractions; tumor-only research profile',fontsize=14,weight='bold');fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(OUT/'regional-common-SNP-allele-balance.png',dpi=170);plt.close(fig)
print(json.dumps([r for r in summ if r['filter']=='candidate_heterozygous_stringent'],indent=2))

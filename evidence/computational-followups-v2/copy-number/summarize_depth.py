"""Descriptive GC/repeat-adjusted off-exon density; no absolute CN calls."""
import pathlib,json,gzip,csv,collections,time
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=pathlib.Path(__file__).resolve().parent
def tsv(p,rows):
 if not rows:return
 with p.open('w') as f:w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
centro={}
with gzip.open(OUT/'public/cytoBandIdeo.txt.gz','rt') as f:
 for ln in f:
  c,s,e,b,st=ln.strip().split('\t')
  if st=='acen':centro.setdefault(c,[]).extend([int(s),int(e)])
centro={c:[min(p),max(p)] for c,p in centro.items()}
def arm(c,s,e):
 a,b=centro.get(c,[0,0])
 if e<=a:return c+'p'
 if s>=b:return c+'q'
 return 'centromere'
def aggregate(kbins):
 rows=[]
 for ch in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']:
  loaded=np.load(OUT/'genome-bins'/f'{ch}.npz');z={k:loaded[k] for k in loaded.files};loaded.close();n=len(z['length'])
  for i in range(0,n,kbins):
   j=min(n,i+kbins);s=i*1000;e=s+int(z['length'][i:j].sum());r={'chrom':ch,'start0':s,'end0':e,'arm':arm(ch,s,e)}
   for key in z:r[key]=int(z[key][i:j].sum())
   rows.append(r)
 return rows
def design(gc,rep):
 g=gc-.4;r=rep-.4
 return np.column_stack([np.ones(len(g)),g,g*g,g*g*g,r,r*r,g*r])
def robustfit(X,y):
 beta=np.linalg.lstsq(X,y,rcond=None)[0]
 for i in range(30):
  res=y-X@beta;scale=max(.1,float(np.median(np.abs(res-np.median(res)))*1.4826));w=np.minimum(1,1.345*scale/np.maximum(np.abs(res),1e-12));b=np.linalg.lstsq(X*np.sqrt(w[:,None]),y*np.sqrt(w),rcond=None)[0]
  if np.max(np.abs(beta-b))<1e-8:break
  beta=b
 return beta
method={'scope':'relative off-exon DNA fragment density, NOT measured copy number','capture_proxy':'GENCODE37 protein-coding exon union with100/500/2000bp padding; not actual capture BED','normalization':'100/250kb bins. Density numerator proper-pair read1 fragment midpoints, denominator off-exon ACGT bp. Robust within-sample regression of log2 density on cubic GC, quadratic repeat fraction and interaction. Train on autosomes excludingchr3/5/9/17, excludingcentromeres. Center to trainingmedian. Genome background is not demonstrated diploid.','filters':'retainedoff-exonACGT>=40%bin;genomeACGT>=90%;GC0.2–0.75;repeat<=0.85;atleast5fragments','pseudocount':'0.5 added to fragment count before log2 to bound zero bins; retainedbins>=5 counts','limitations':['Same-sample GC/RepeatMasker correction cannot remove unknown bait/mappability/protocol biases.','No matched normal/panel of normals/assay target BED/purity/ploidy; relative density is not absolute or allele-specific copy number.','Off-exon filter uses fragment center, not full fragment; flank sensitivity included.','Centromere and sex chromosomes excluded from normalization.','Nearby windows and markers are correlated; plotting/region medians do not create independent molecular observations.'],'fits':{}}
allrows=[];summaries=[]
for kbins in [100,250]:
 rows=aggregate(kbins)
 for pad,mq in [(100,30),(500,30),(2000,30),(500,60)]:
  tag=f'off{pad}'+('_mq60' if mq==60 else '');den=np.array([r[f'off{pad}_acgt_bases'] for r in rows],float);gc=np.array([r[f'off{pad}_gc_bases'] for r in rows])/np.maximum(den,1);rep=np.array([r[f'off{pad}_repeat_bases'] for r in rows])/np.maximum(den,1);count=np.array([r['fragment_midpoints_'+tag] for r in rows],float);length=np.array([r['length'] for r in rows],float);acgt=np.array([r['acgt_bases'] for r in rows]);logdens=np.log2((count+.5)/np.maximum(den,1)*1000)
  keep=(den>=.4*length)&(acgt>=.9*length)&(gc>=.2)&(gc<=.75)&(rep<=.85)&(count>=5)&np.array([r['arm']!='centromere' for r in rows]);train=keep&np.array([r['chrom'] not in ['chr3','chr5','chr9','chr17','chrX','chrY'] for r in rows]);X=design(gc,rep);beta=robustfit(X[train],logdens[train]);res=logdens-X@beta;center=float(np.median(res[train]));res-=center;rawcenter=float(np.median(logdens[train]));fitkey=f'{kbins}kb_{tag}';method['fits'][fitkey]={'training_bins':int(train.sum()),'retained_bins':int(keep.sum()),'beta':beta.tolist(),'residual_center':center,'raw_log2density_center':rawcenter}
  for i,r in enumerate(rows):
   r[f'{tag}_usable']=int(keep[i]);r[f'{tag}_gc']=float(gc[i]);r[f'{tag}_repeat_fraction']=float(rep[i]);r[f'{tag}_raw_log2_relative_density']=float(logdens[i]-rawcenter);r[f'{tag}_adjusted_log2_relative_density']=float(res[i])
  for a in sorted(set(r['arm'] for r in rows)):
   ix=np.array([r['arm']==a for r in rows])&keep
   if not ix.any():continue
   summaries.append({'bin_kb':kbins,'mask_flank':pad,'MAPQ_min':mq,'arm':a,'bins':int(ix.sum()),'raw_log2_density_median':float(np.median(logdens[ix]-rawcenter)),'adjusted_log2_density_median':float(np.median(res[ix])),'adjusted_density_ratio_median':float(2**np.median(res[ix])),'adjusted_log2_Q25':float(np.quantile(res[ix],.25)),'adjusted_log2_Q75':float(np.quantile(res[ix],.75)),'retained_off_exon_fragments':int(count[ix].sum()),'retained_off_exon_ACGT_bases':int(den[ix].sum())})
 tsv(OUT/f'genome-relative-density.{kbins}kb.tsv',rows);allrows.extend(rows if kbins==100 else [])
tsv(OUT/'chromosome-arm-density-summary.tsv',summaries)
(OUT/'relative-depth-normalization-method.json').write_text(json.dumps(method,indent=2)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'font.size':10})
fig,axes=plt.subplots(4,1,figsize=(12,11),sharey=True)
for ax,ch in zip(axes,['chr3','chr5','chr9','chr17']):
 r=[x for x in allrows if x['chrom']==ch];xx=np.array([(x['start0']+x['end0'])/2e6 for x in r]);valid=np.array([x['off500_usable'] for x in r],bool);yy=np.array([x['off500_adjusted_log2_relative_density'] for x in r]);ax.scatter(xx[valid],yy[valid],s=4,color='#888888',alpha=.45,label='100kb bins, exon+500bp exclusion')
 for pad,color in [(100,'#111111'),(500,'#c34b36'),(2000,'#28747e')]:
  med=[];xs=[]
  for start in range(0,int(max(xx))+1,2):
   vv=[x[f'off{pad}_adjusted_log2_relative_density'] for x in r if start<=x['start0']/1e6<start+2 and x[f'off{pad}_usable']]
   if len(vv)>=4:
    if xs and start+1-xs[-1]>2.1:xs.append((start+1+xs[-1])/2);med.append(np.nan)
    xs.append(start+1);med.append(np.median(vv))
  ax.plot(xs,med,lw=1.4,color=color,label=f'2Mb median, exon+{pad}bp exclusion')
 ax.axhline(0,color='#444444',ls=':',lw=.8);a,b=centro[ch];ax.axvspan(a/1e6,b/1e6,color='#dddddd',alpha=.6);ax.set_title(ch,loc='left',weight='bold');ax.set_ylim(-2.2,2.2);ax.set_ylabel('Relative log₂ density');ax.set_xlabel('Genomic position (Mb, hg38)')
 for gene,pos in [('BAP1',52.40855),('RASA1',87.33256),('APC',112.843324),('MTAP / CDKN2A/B',21.9)]:
  if (ch=='chr3' and gene=='BAP1') or (ch=='chr5' and gene in ['RASA1','APC']) or (ch=='chr9' and gene.startswith('MTAP')):ax.axvline(pos,color='#111111',ls='--',lw=.9);ax.text(pos+.8,1.85,gene,fontsize=9)
axes[0].legend(loc='lower left',ncol=2,fontsize=8);fig.suptitle('DNA outside annotated exons: relative density and mask sensitivity\nSingle-sample research profile — not calibrated copy number',fontsize=14,weight='bold');fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(OUT/'regional-off-exon-relative-density.png',dpi=170);plt.close(fig)
print(json.dumps([r for r in summaries if r['bin_kb']==100 and r['mask_flank']==500 and r['MAPQ_min']==30 and r['arm'].startswith(('chr3','chr5','chr9','chr17'))],indent=2))

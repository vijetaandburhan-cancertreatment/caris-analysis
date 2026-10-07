from pathlib import Path
import csv,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
rows=[r for r in csv.DictReader((OUT/'candidate-gtex-comparison.tsv').open(),delimiter='\t') if r['reference_aware_research_pass']=='1' and not r['GTEx_samples_with_junction']]
evidence=json.loads((OUT/'candidate-reference-read-evidence.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
fig,axs=plt.subplots(len(rows),1,figsize=(10,8.5));fig.subplots_adjust(left=.16,right=.96,top=.88,bottom=.11,hspace=.65)
for ax,r in zip(axs,rows):
 key=f"{r['gene']}:{r['chrom']}:{r['intron_start0']}-{r['intron_end0']}";ev=evidence[key];uniq={}
 for e in ev['examples']:
  if e['perfect12'] and (e['query_name'] not in uniq or min(e['left_anchor'],e['right_anchor'])>min(uniq[e['query_name']]['left_anchor'],uniq[e['query_name']]['right_anchor'])):uniq[e['query_name']]=e
 ex=list(uniq.values())[:9]
 for y,e in enumerate(ex):
  ax.plot([-min(90,e['left_anchor']),0],[y,y],color='black',lw=3);ax.plot([65,65+min(90,e['right_anchor'])],[y,y],color='black',lw=3);ax.plot([0,65],[y,y],color='0.5',ls='--',lw=.7)
 ax.axvline(0,color='0.6',lw=.8);ax.axvline(65,color='0.6',lw=.8)
 ax.set_xlim(-95,160);ax.set_ylim(-1,max(3,len(ex)))
 ax.set_yticks(range(len(ex)));ax.set_yticklabels([f"Example {i+1}" for i in range(len(ex))],fontsize=8)
 ax.set_xticks([-80,-40,0,65,105,145]);ax.set_xticklabels(['−80','−40',r['intron_start0'],r['intron_end0'],'+40','+80'],fontsize=8)
 ax.set_title(f"{r['gene']}  {r['chrom']}  |  {r['fragments']} fragment names, motif {r['reference_splice_motif']}\n{float(r['fraction_same_left_boundary'])*100:.2f}% of same-left-boundary junction fragments; no exact GTEx match",loc='left',weight='bold',fontsize=11)
 ax.set_xlabel('Genomic coordinates increase left to right; both genes are on the minus strand. Intron compressed.',fontsize=8)
 for k in ['top','right']:ax.spines[k].set_visible(False)
 ax.text(.5,.98,'Compressed intron',transform=ax.transAxes,ha='center',va='top',fontsize=8,backgroundcolor='white')
fig.suptitle('Two low-level RNA junctions retained for research review',fontsize=15,weight='bold',y=.97)
fig.text(.08,.025,'Aligned anchor lengths from the existing BAM; one example per query name. Full CIGARs and anchor bases are saved.\nAbsence from the queried GTEx atlas does not establish tumor specificity. Nearby DNA reads were predominantly reference.\nNo clinical splice-mutation call, drug-response prediction or tumor-specific antigen claim was made.',fontsize=8)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'rare-junction-read-evidence.{ext}',dpi=180)
print(OUT/'rare-junction-read-evidence.png')

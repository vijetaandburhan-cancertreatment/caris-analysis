from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,PathPatch
from matplotlib.path import Path as MPath
OUT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.dpi':150,'axes.titlesize':13})
fig,axs=plt.subplots(3,1,figsize=(9.5,10.5));fig.subplots_adjust(hspace=.35,top=.92,bottom=.075,left=.06,right=.98)
def setup(ax,title):
 ax.set_xlim(-.8,6.6);ax.set_ylim(-1.1,2.1);ax.axis('off');ax.set_title(title,loc='left',weight='bold',pad=12)
def exon(ax,x,text,face='0.9',w=1):
 ax.add_patch(Rectangle((x,-.1),w,.25,facecolor=face,edgecolor='black',lw=1));ax.text(x+w/2,-.4,text,ha='center',va='top',fontsize=10)
def arc(ax,x1,x2,y,label,dashed=False):
 path=MPath([(x1,.15),(x1,y),(x2,y),(x2,.15)],[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4])
 ax.add_patch(PathPatch(path,fill=False,lw=1.3,ls='--' if dashed else '-',color='black'))
 ax.text((x1+x2)/2,y*.74+.15,label,ha='center',va='bottom',fontsize=10,backgroundcolor='white')
setup(axs[0],'A. MET: no observed exon-14-skipping junction')
for x,t in [(0,'Exon 13'),(2.3,'Exon 14'),(4.6,'Exon 15')]:exon(axs[0],x,t)
arc(axs[0],1,2.3,.65,'88 fragments');arc(axs[0],3.3,4.6,.65,'87 fragments');arc(axs[0],1,4.6,1.65,'Exon 13 → 15: 0 reads',True)
axs[0].text(-.5,-.86,'Zero skipping support also holds before anchor filtering; inclusion has 169 / 150 reads.\nThis is a bounded negative for this exact junction, not a validated clinical exclusion.',va='top',fontsize=9)
setup(axs[1],'B. BAP1: the frameshift is physically linked to canonical splicing')
for x,t in [(0,'Exon 3'),(2.3,'Exon 4\ncontains c.168_178del11'),(4.6,'Exon 5')]:exon(axs[1],x,t)
arc(axs[1],1,2.3,.72,'141 direct-read links');arc(axs[1],3.3,4.6,.72,'31 direct-read links')
arc(axs[1],1,4.6,1.8,'Exon-4 skip: 2 fragments; junction in an annotated NMD isoform',True)
axs[1].plot([2.8,2.8],[-.11,.17],color='black',lw=3)
axs[1].text(-.5,-.89,'Direct links require the deletion and junction in the same read (mates collapsed by name).\nLocal mutant-spliced RNA support does not establish a full isoform, protein or HLA display.',va='top',fontsize=9)
setup(axs[2],'C. APC: the strongest GENCODE37-unannotated event is a natural isoform')
exon(axs[2],0,'Upstream exon');exon(axs[2],3,'Full exon 10*',w=2.3)
axs[2].add_patch(Rectangle((3,-.1),1.5,.25,facecolor='0.5',edgecolor='black',hatch='///'))
arc(axs[2],1,3,.62,'226 fragments');arc(axs[2],1,4.5,1.75,'222 fragments: 303-nt partial-exon exclusion')
axs[2].text(-.5,-.84,'This matches the established APC exon-9a splice form, r.934_1236del.\n*Exon 10 in the selected transcript includes the historical coding-exon-9 region.',va='top',fontsize=9)
fig.suptitle('RNA splice-junction evidence: useful results and false-lead checks',fontsize=15,weight='bold',y=.975)
fig.text(.06,.022,'Schematic spacing; exon numbering follows the selected MANE transcripts. Fragment counts are read-group + paired-read names, not UMIs.\nMain inventory: primary, MAPQ≥20, NH=1 where present; ≥12-base BQ20 anchors on both sides. Source: Caris RNA BAM, GENCODE37.',fontsize=8)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'key-splice-evidence.{ext}',dpi=180)
print(OUT/'key-splice-evidence.png')

from pathlib import Path
import json,csv,html,sys
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
D=Path(__file__).resolve().parent;R=D.parent;preview='--preview' in sys.argv
if not preview:assert json.loads((D/'final-content-approval.json').read_text())['approved']
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
old=rows(R/'frozen-exception-sites.tsv');new={(x['chrom'],x['pos1']):x for x in rows(R/'comparison/all16-exception-comparison.tsv')}
s={k:ParagraphStyle(k,fontName=f,fontSize=z,leading=l,spaceAfter=a) for k,f,z,l,a in [('body','Times-Roman',10.2,13,6),('small','Helvetica',8.3,10.5,5),('h','Helvetica-Bold',11,14,6),('title','Helvetica-Bold',18,22,7),('cell','Times-Roman',9.1,11.6,0),('th','Helvetica-Bold',8.8,11.2,0)]}
s['h'].keepWithNext=True;s['h'].spaceBefore=6
story=[];plain=[]
def para(t,style='body'):
 plain.append(t);return Paragraph(html.escape(t).replace('\n','<br/>'),s[style])
def add(t,style='body'):story.append(para(t,style))
def table(rs,widths):
 t=Table([[para(str(v),'th' if i==0 else 'cell') for v in row] for i,row in enumerate(rs)],colWidths=widths,hAlign='LEFT',repeatRows=1);t.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('LINEABOVE',(0,0),(-1,0),.7,colors.black),('LINEBELOW',(0,0),(-1,0),.5,colors.black),('LINEBELOW',(0,-1),(-1,-1),.5,colors.black),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]));story.extend([t,Spacer(1,6)])
def ac(vals):return ', '.join(b+str(n) for b,n in zip('ACGT',vals) if int(n)) or '0'
add('Full-pair RNA alignment comparison','title');add('Caris TN26-279853 | Fixed-cohort research QC | 5 October 2026','small')
add('Five of the six stronger DNA/RNA discrepancy positions persist under the second alignment pipeline, including both linked PABPC1 sites. The APOD discrepancy is absent in the new counts and is demonstrably alignment-dependent. The experiment does not establish which placement is biologically correct or change the existing treatment-related research priorities.')
add('Same starting cohort; coverage losses retained','h')
add('All 23,209,264 original RNA read pairs were processed against the resident full primary-genome reference. The comparison was frozen at 4,015 originally jointly callable public SNPs, with the original DNA counts and categories unchanged. The new RNA remains callable at depth ≥20 at 3,932 sites; 83 fall below 20, and none loses all A/C/G/T evidence.')
table([['Among the same 2,873 DNA-dominant sites','Original RNA','New RNA'],['Depth ≥20 and same allele ≥90%','2,860 / 2,873','2,805 / 2,873'],['Depth ≥20 and same allele ≥98%','2,760 / 2,873','2,721 / 2,873'],['Depth below 20 (kept in denominator)','0','59'],['Depth ≥20 but same allele <90%','13','9']],[315,104,104])
add('Of the original 13 exceptions, five reach ≥90% same-allele retention and eight remain below it. One additional, low-depth MLLT3 site crosses below 90%. The 59 coverage failures came from previously passing sites. No baseline DNA-dominant site has a callable ≥98% opposite REF/ALT allele. Fraction-only retention that ignores minimum depth is a separate statistic and is not substituted for the table above.','small')
add('Six stronger original discrepancies','h')
rs=[['Locus (GRCh38, 1-based)','Original DNA','Original RNA','New RNA']]
for o in old:
 if o['review_class']!='strong_baseline_screen':continue
 n=new[(o['chrom'],o['pos1'])];rs.append([o['overlap_gene_symbols'].replace(';',' / ')+'\n'+o['chrom']+':'+o['pos1'],ac(o['baseline_DNA_ACGT'].split(',')),ac(o['baseline_RNA_ACGT'].split(',')),ac([n['new_'+b] for b in 'ACGT'])])
table(rs,[184,80,132,127])
add('Counts are unambiguous query-name units, not unique molecules. The strong screen requires both DNA/RNA depths ≥50 and a base at ≤2% in DNA but ≥10% and at least 10 units in RNA. APOD no longer meets it; the other five still do.','small')
add('PABPC1 remains unresolved','h')
add('The new data retain G at chr8:100705591 in 61/470 names and C at chr8:100705604 in 73/430. All 61 G-supporting names also support C on the same qualifying aligned records; 56 of the original 65 shared names still support both. These linked observations must not be treated as two independent checks or converted into tumor mutations. Their origin remains unresolved.')
story.append(PageBreak())
add('Read context and remaining exceptions','title')
add('Complete emitted records were preserved for all 3,426 preselected discrepancy-associated query names, including mates, secondary/supplementary placements and unmapped records. All names were recovered. New capture retained 9,061 named records, not merely alignments overlapping the questioned positions.')
add('For APOD, 26 of 27 original C-supporting names now have qualifying G evidence at the site; the remaining name overlaps it but fails matched filters. For EXOC3, 14 of 22 original A names retain A; three now support C, one conflicts, two fail filters and two are unmapped. For MSH3/DHFR, 58 of 104 original G names retain G; 41 map without a base at this position, three fail filters and two are unmapped. These are changed alignment/filter outcomes, not proof of corrected origin.')
add('Other original exceptions remain in the record','h')
rs=[['Group / locus','Original RNA','New RNA']]
for o in old:
 if o['review_class']=='strong_baseline_screen':continue
 n=new[(o['chrom'],o['pos1'])];tag='MAPQ60-only: ' if o['review_class']=='MAPQ60_only_not_baseline_flag' else ''
 rs.append([tag+o['overlap_gene_symbols'].replace(';',' / '),ac(o['baseline_RNA_ACGT'].split(',')),ac([n['new_'+b] for b in 'ACGT'])])
table(rs,[262,132,129])
add('The three MAPQ60-only flags were originally produced by tightening the DNA mapping-quality filter; they remain separate sensitivity findings. This table shows their matched baseline-filter RNA counts. MLLT3 chr9:20414378 is newly below 90%: original A24/G2 becomes A24/G3 (depth 27), below the strong screen’s depth requirement.','small')
add('Methods and verification','h')
add('Matched count filters require MAPQ ≥30, BQ ≥25, proper pairing and exclusion of flags 0xF0C. Qualifying A/C/G/T observations are collapsed by exact query name; conflicting bases are excluded. Coverage means qualifying A/C/G/T, not reference skips. Capture preserves all SNP-base-overlapping records plus every emitted record for frozen names; D/N-only bookkeeping for other names is incomplete by design.')
add('An independent aligned-pair implementation reproduced all 4,015 per-base count rows and material exception transitions. The capture logic passed synthetic boundary/flag/name controls; final BAM quick-check and complete input-pair counts passed. The complete record and per-name tables, scripts, frozen design and audit receipts are retained in the evidence archive.')
add('Original Caris STAR 2.7.8a/CTAT/two-pass settings differ from local STAR 2.7.11b/primary-genome/D8/one-pass settings, including overlap, splice and chimeric parameters. This is not a controlled sparsity-only experiment, a validated identity/contamination assay or a clinical sensitivity test. New placements are not automatically corrections. The clean-read sensitivity and all 16 original review sites remain available, including coverage losses.')
add('Local full pass completed in 26.8 minutes. The selective BAM is 64.1 MB; monitored peak main-pipeline RSS was 8.40 GB and minimum free disk 4.99 GB. No full BAM, new clinical target or external outreach resulted. The biological v2 handoff remains unchanged.','small')
def foot(c,d):
 c.saveState();w,h=A4;c.setStrokeColor(colors.black);c.setLineWidth(.4);c.line(36,29,w-36,29);c.setFont('Helvetica',7.5);c.drawString(36,18,'TN26-279853 | Alignment-method comparison | 5 October 2026');c.drawRightString(w-36,18,str(d.page));c.restoreState()
name='Caris-final-alignment-comparison-2026-10-05'+('-preview' if preview else '');out=D/(name+'.pdf')
SimpleDocTemplate(str(out),pagesize=A4,leftMargin=36,rightMargin=36,topMargin=33,bottomMargin=39,title='Caris final alignment-method comparison',author='Research coordination').build(story,onFirstPage=foot,onLaterPages=foot)
(D/(name+'.txt')).write_text('\n\n'.join(plain)+'\n');print(out)

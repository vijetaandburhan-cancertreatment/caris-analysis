from pathlib import Path
import json,csv,html,hashlib,sys
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
D=Path(__file__).resolve().parent;R=D.parent
x=json.loads((D/'editorial.json').read_text())
preview='--preview' in sys.argv
assert preview or x['ready_for_final_render'] is True
s={k:ParagraphStyle(k,fontName=font,fontSize=size,leading=lead,spaceAfter=after) for k,font,size,lead,after in [('body','Times-Roman',10.7,13.5,7),('small','Helvetica',8.4,10.7,6),('h','Helvetica-Bold',11.4,14,6),('title','Helvetica-Bold',19,23,7),('cell','Times-Roman',9.1,11.6,0),('th','Helvetica-Bold',8.9,11.3,0)]}
s['h'].keepWithNext=True;s['h'].spaceBefore=6
story=[];plain=[]
def p(t,sty='body'):
 plain.append(t);return Paragraph(html.escape(t).replace('\n','<br/>'),s[sty])
def add(t,sty='body'):story.append(p(t,sty))
def tab(rows,widths):
 z=[[p(c,'th' if i==0 else 'cell') for c in row] for i,row in enumerate(rows)]
 t=Table(z,colWidths=widths,hAlign='LEFT',repeatRows=1);t.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('LINEABOVE',(0,0),(-1,0),.7,colors.black),('LINEBELOW',(0,0),(-1,0),.5,colors.black),('LINEBELOW',(0,-1),(-1,-1),.5,colors.black),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LEFTPADDING',(0,0),(-1,-1),5)]));story.extend([t,Spacer(1,7)])
add('DNA/RNA consistency check','title');add('Caris TN26-279853 | Population-selected SNP panel | Technical addendum | 5 October 2026','small')
add('Broad allele compatibility was reproduced across a larger panel selected without using this patient’s DNA allele calls. Local discrepancies remain. These results support continued research use of the paired data; they do not establish sample identity, quantify contamination, validate a treatment target, or change the prior biological findings.')
add('What was tested','h')
add('The frozen public panel contains 165,782 common SNPs on chromosomes 2, 3, 5, 8, 9 and 17. Selection used population frequencies and preselected exon, regional and genomic-grid intervals. Every position remains in the accounting, including those with no qualifying bases. A callable position requires at least 20 unambiguous A/C/G/T query-name units in each assay.')
tab([['Coverage stage','DNA','RNA'],['At least one qualifying A/C/G/T unit','62,262','25,707'],['At least 20 qualifying A/C/G/T units','35,885','4,072'],['At least 20 in both assays','4,015 shared positions','4,015 shared positions']],[295,114,114])
add('Result at positions dominated by one DNA allele','h')
add('Of 4,015 jointly callable positions, 2,873 have one REF or ALT allele accounting for at least 98% of DNA units. RNA carries that same allele at ≥90% at 2,860/2,873 positions (99.55%), and at ≥98% at 2,760/2,873 (96.07%). None of these 2,873 baseline DNA-dominant positions switches to ≥98% of the opposite REF/ALT allele.')
tab([['DNA category','Positions','RNA same allele ≥90%','RNA same allele ≥98%'],['REF dominant','2,283','2,272','2,189'],['ALT dominant','590','588','571'],['Total','2,873','2,860','2,760']],[181,76,133,133])
add('Breadth and selection limits','h')
add('The dominant positions span all six chromosomes, 554 one-megabase bins and 1,194 uniquely assigned GENCODE genes. Only 72 overlap the earlier 335-marker panel. They are still clustered: 549 adjacent gaps are <1 kb, and surviving sites are strongly enriched near annotated exons. They must not be treated as independent random observations.')
add('A coordinate-first subset spaced ≥50 kb apart was frozen independently of allele direction. Of its 17,264 selected sites, 165 were jointly callable and 120 were DNA dominant: all 120 retained the same RNA allele at ≥90%; 115 did so at ≥98%. This supports spatial breadth, not a probability of identity.')
add('The remaining 1,142 jointly callable positions comprise 1,109 mixed REF/ALT DNA sites and 33 residual mixtures. RNA becomes ≥98% REF or ALT at 68 of the mixed-DNA sites. Such allele-fraction shifts are reported separately and are not counted as identity mismatches.','small')
story.append(PageBreak())
add('Discrepancies retained for review','title')
add('The prespecified strong-discrepancy screen required ≥50 A/C/G/T units in each assay, an individual base at ≤2% in DNA but ≥10% and at least 10 units in RNA. Independent per-base counting across all 4,015 sites confirmed six positions. Counts below are query-name units, not unique molecules. Coordinates are GRCh38, 1-based.','body')
tab([['Locus','DNA baseline','RNA baseline','RNA after clean-read filter'],['COL5A2\nchr2:189043211','G58','G327, A116','G183, A0'],['APOD\nchr3:195568834','G87','G36, C27','G0, C24'],['EXOC3\nchr5:443259','C153','C74, A22, G1','C11, A21, G1'],['MSH3 / DHFR\nchr5:80654917','C912','C9, G104','C7, G53'],['PABPC1\nchr8:100705591','A94','A400, G65','A278, G9'],['PABPC1\nchr8:100705604','G101','G358, C78, A1','G218, C31']],[176,67,136,144])
add('The clean-read sensitivity excludes soft-clipped reads and requires at least five query bases between the tested base and each aligned end. It changes both evidence and coverage. COL5A2’s unexpected A support disappears; APOD and EXOC3 remain discrepant but fall below depth 50. Falling below a threshold is not resolution. The two PABPC1 positions are only 13 bases apart.','small')
add('Read-context review','h')
for t in x['origin_paragraphs']:add(t)
add('Other exceptions and mapping sensitivity','h')
add('All 13 dominant-DNA sites with <90% same-allele RNA retention are preserved in the evidence table, including seven below the strong screen’s depth requirement. Raising MAPQ to 60 introduces three further flags—FANCD2, SDHAP2/MUC20-OT1 and FOXD4L4—by removing DNA allele support while leaving RNA unchanged. These are sensitivity-only findings, not replacements for the baseline results.')
add('No discrepancy has been converted into a new tumor mutation, neoantigen, diagnosis or drug recommendation. Locus-specific biology and mapping effects are not separated by the aggregate consistency result.','small')
story.append(PageBreak())
add('Methods, verification and scope','title')
add('Frozen before counting','h')
add('Design SHA256: 90a153f630375103363792037618d93a864b01bad029575f36d6a3a382ee0bfa','small')
add('Population selection used 1000 Genomes 2019 GRCh38 biallelic PASS variants (global ALT frequency 0.05–0.95), protein-coding gene exons plus 100 bp, designated gene regions plus 250 kb, and a common-SNP grid. The six-chromosome scope was inherited from the prior copy-number work. SAS frequency was recorded, not used for selection. Every selected REF matched the resident GRCh38 FASTA.')
add('Both original Caris BAMs were recounted with samtools 1.24 / pysam 0.24.1. Filters: MAPQ ≥30; base quality ≥25; proper pairs required; flags 0xF0C excluded; BAQ and samtools mate-overlap adjustment disabled. Qualifying observations were collapsed by query name; conflicting A/C/G/T mates were discarded. Each original BAM has one read group. Reference skips, deletions and ambiguous bases do not count toward A/C/G/T depth.')
add('The explicit depth limit was disabled (-d 0, effectively unlimited in samtools). Installed behavior was qualified with 9,001 synthetic paired names / 18,002 records, a finite-cap contrast, and separate base-quality, mapping-quality, duplicate, supplementary and conflicting-mate controls. Pileup columns and BED coordinate conventions were checked before patient counting.')
add('Independent checks','h')
add('A separate CIGAR-walk implementation reproduced all 4,015 jointly callable DNA and RNA positions: REF/ALT/other query-name counts, conflicts and REF/ALT read counts. Exact integer arithmetic independently reproduced categories, retention thresholds, full-panel attrition and coordinate-first spacing. Public GENCODE annotation passed an independent 81-position boundary/category check.')
add('The fresh uncapped DNA run exactly matched all 935,370 compared integer fields across the 62,358 previously observed DNA rows. Thus the old finite cap had no observed effect on these final counts. This conclusion comes from an actual uncapped recount, not the fact that reported post-quality depth was below the old cap.')
add('Newer RNA alignment: not evaluable here','h')
add('The newer full-reference STAR run was streamed and did not retain a complete BAM. Its retained BAP1/RASA1/LATS1/LATS2 BAM cannot substitute for the full panel. Only 73 selected public SNPs intersect BAP1/RASA1 retained intervals; none met joint callability (maximum original RNA depth 9). No per-site newer-alignment comparison was performed, and thresholds were not relaxed. The broad result therefore depends on the original Caris alignments.')
add('Interpretation and files','h')
add('This is a research consistency analysis on selected, coverage-surviving loci—not a validated identity assay, chance-match calculation, contamination estimate, matched-normal analysis, or proof of allele origin. Query names are not UMIs. The original four follow-ups and v2 handoff remain unchanged; this addendum adds technical quality evidence, not a treatment recommendation.')
add('The accompanying evidence archive includes the frozen design, counts for all selected loci, coverage and discrepancy tables, independent audit code/results, read-context evidence, source provenance and hashes. Large compressed raw pileups remain in the resident analysis folder. Local execution only; no new AWS compute resources or external outreach.','small')
add('Method reference: https://www.htslib.org/doc/samtools-mpileup.html (samtools 1.24). GENCODE37 annotation and 1000 Genomes panel provenance are recorded in the evidence files.','small')
def foot(c,d):
 c.saveState();w,h=A4;c.setStrokeColor(colors.black);c.setLineWidth(.4);c.line(36,29,w-36,29);c.setFillColor(colors.black);c.setFont('Helvetica',7.5);c.drawString(36,18,'TN26-279853 | Research quality-control addendum | 5 October 2026');c.drawRightString(w-36,18,str(d.page));c.restoreState()
out=D/('Caris-DNA-RNA-population-panel-addendum-2026-10-05'+('-preview' if preview else '')+'.pdf')
SimpleDocTemplate(str(out),pagesize=A4,leftMargin=36,rightMargin=36,topMargin=33,bottomMargin=39,title='Caris DNA/RNA population-panel consistency addendum',author='Research coordination').build(story,onFirstPage=foot,onLaterPages=foot)
(D/(out.stem+'.txt')).write_text('\n\n'.join(plain)+'\n');print(out)

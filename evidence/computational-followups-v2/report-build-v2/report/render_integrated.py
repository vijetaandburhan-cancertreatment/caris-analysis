from pathlib import Path
import json,html,shutil,hashlib
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
W=Path(__file__).resolve().parent;O=Path('outputs/caris-followup');D=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct5-followup-findings-v2');D.mkdir(parents=True,exist_ok=True)
x=json.loads((W/'report-content.json').read_text());fusion=json.loads((W/'fusion-status.json').read_text());execution=json.loads((W/'execution-status.json').read_text());stem='Caris-followup-findings-2026-10-05-v2';out=O/(stem+'.pdf')
s={'body':ParagraphStyle('body',fontName='Times-Roman',fontSize=10.7,leading=13.5,spaceAfter=7),'h':ParagraphStyle('h',fontName='Helvetica-Bold',fontSize=11.4,leading=14,spaceBefore=7,spaceAfter=6,keepWithNext=True),'title':ParagraphStyle('title',fontName='Helvetica-Bold',fontSize=19,leading=23,spaceAfter=5),'page':ParagraphStyle('page',fontName='Helvetica-Bold',fontSize=15,leading=18,spaceAfter=9),'small':ParagraphStyle('small',fontName='Helvetica',fontSize=8.5,leading=11,spaceAfter=7),'cell':ParagraphStyle('cell',fontName='Times-Roman',fontSize=10,leading=12,spaceAfter=0),'th':ParagraphStyle('th',fontName='Helvetica-Bold',fontSize=9.3,leading=12)}
def p(t,sty='body'):return Paragraph(html.escape(t).replace('\n','<br/>'),s[sty])
story=[];plain=[x['title'],x['subtitle'],x['status']]
for i,page in enumerate(x['pages']):
 if i:story.append(PageBreak())
 if i==0:story.extend([p(x['title'],'title'),p(x['subtitle'],'small'),p(x['status'],'small')])
 story.extend([p(page['title'],'page'),p(page['intro'])]);plain.extend([page['title'],page['intro']])
 for sec in page['sections']:
  story.append(p(sec['heading'],'h'));plain.append(sec['heading'])
  if 'text' in sec:story.append(p(sec['text']));plain.append(sec['text'])
  if 'table' in sec:
   tab=sec['table'];widths=[150,81,81] if len(tab[0])==3 else [115,408]
   if len(tab[0])==3:widths=[290,116,117]
   data=[[p(c,'th' if r==0 else 'cell') for c in row] for r,row in enumerate(tab)]
   table=Table(data,colWidths=widths,hAlign='LEFT',repeatRows=1);table.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('LINEABOVE',(0,0),(-1,0),.7,colors.black),('LINEBELOW',(0,0),(-1,0),.5,colors.black),('LINEBELOW',(0,-1),(-1,-1),.5,colors.black),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LEFTPADDING',(0,0),(-1,-1),5)]));story.extend([table,Spacer(1,7)]);plain.append('\n'.join(' | '.join(r) for r in tab))
  if 'after' in sec:story.append(p(sec['after']));plain.append(sec['after'])
 if i==len(x['pages'])-1:
  story.extend([p(fusion['heading'],'h'),p(fusion['text'])]);plain.extend([fusion['heading'],fusion['text']])
  story.extend([p(execution['heading'],'h'),p(execution['text'],'small')]);plain.extend([execution['heading'],execution['text']])
  sources='Evidence: copy-number; hla-II; HLA allele support; validation package; targeted splicing; MET exon-14; fusion read-level audit, public robustness controls, CCND1::HEG1 RNA/DNA audits and diagnostic counters. Methods: github.com/changyunjian/CapHLA; github.com/suhrig/arriba; doi.org/10.1093/bib/bbae595. Research only; originals preserved.'
  story.append(p(sources,'small'));plain.append(sources)
def foot(c,d):
 c.saveState();w,h=A4;c.setLineWidth(.4);c.line(36,29,w-36,29);c.setFont('Helvetica',7.6);c.drawString(36,18,'Computational findings | Research interpretation | 5 October 2026');c.drawRightString(w-36,18,str(d.page));c.restoreState()
doc=SimpleDocTemplate(str(out),pagesize=A4,leftMargin=36,rightMargin=36,topMargin=33,bottomMargin=39,title=x['title'],author='Research coordination')
doc.build(story,onFirstPage=foot,onLaterPages=foot)
(O/(stem+'.txt')).write_text('\n\n'.join(plain)+'\n')
for q in [out,O/(stem+'.txt')]:shutil.copy2(q,D/q.name)
for q in W.glob('*.json'):shutil.copy2(q,D/q.name)
shutil.copy2(W/'render_integrated.py',D/'render_integrated.py');
if (W/'audit_summary.py').exists():shutil.copy2(W/'audit_summary.py',D/'audit_summary.py')
manifest=[]
for q in sorted(D.iterdir()):
 if q.is_file() and q.name!='checksums.sha256':manifest.append(hashlib.sha256(q.read_bytes()).hexdigest()+'  '+q.name)
(D/'checksums.sha256').write_text('\n'.join(manifest)+'\n')
print(out);print('Mirrored',D);print('Words',len(' '.join(plain).split()))

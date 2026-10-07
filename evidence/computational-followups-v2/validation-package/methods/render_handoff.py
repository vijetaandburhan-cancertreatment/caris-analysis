from pathlib import Path
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
import re,html
O=Path('outputs/caris-followup/validation-package');P=O/'BAP1-LATS-validation-handoff.pdf'
body=ParagraphStyle('body',fontName='Times-Roman',fontSize=10.7,leading=13.4,spaceAfter=7,textColor=colors.black)
head=ParagraphStyle('heading',fontName='Helvetica-Bold',fontSize=12.6,leading=16,spaceBefore=9,spaceAfter=9,textColor=colors.black,keepWithNext=True)
title=ParagraphStyle('title',fontName='Helvetica-Bold',fontSize=21,leading=25,spaceAfter=7)
sub=ParagraphStyle('sub',fontName='Helvetica',fontSize=9.3,leading=13,spaceAfter=17)
small=ParagraphStyle('small',parent=body,fontSize=9.2,leading=12,spaceAfter=5)
parts=(O/'validation-handoff.txt').read_text().split('\n\n');story=[]
for i,b in enumerate(parts):
 lines=b.splitlines()
 if i==0:
  story.append(Paragraph(html.escape(lines[0]),title));story.append(Paragraph(html.escape(lines[1]),sub));continue
 if re.match(r'^\d\. ',b) or b in ('Package guide','Sources and provenance'):
  story.append(Paragraph(html.escape(b),head));continue
 style=small if b.startswith('[1]') or b.startswith('variants.tsv:') else body
 text=html.escape(b).replace('\n','<br/>')
 text=re.sub(r'(https://[^ &<;]+)',lambda m:f'<link href="{m[1]}" color="black">{m[1]}</link>',text)
 story.append(Paragraph(text,style))
def footer(c,d):
 c.saveState();w,h=A4;c.setStrokeColor(colors.black);c.setLineWidth(.4);c.line(43,38,w-43,38);c.setFont('Helvetica',8);c.drawString(43,25,'BAP1 / LATS research validation | 4 October 2026');c.drawRightString(w-43,25,str(d.page));c.restoreState()
doc=SimpleDocTemplate(str(P),pagesize=A4,rightMargin=43,leftMargin=43,topMargin=38,bottomMargin=49,title='BAP1 and LATS research validation handoff',author='Research coordination',pageCompression=1)
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(P)

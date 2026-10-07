from pathlib import Path
from reportlab.platypus import SimpleDocTemplate,Paragraph
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader
import re,json,shutil,hashlib
W=Path(__file__).resolve().parent;ROOT=W.parents[2];O=ROOT/'output/pdf';O.mkdir(exist_ok=True,parents=True);D=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct5-raw-caris-handoff-v2');D.mkdir(parents=True,exist_ok=True)
pdfmetrics.registerFont(TTFont('Arial','/System/Library/Fonts/Supplemental/Arial.ttf'));pdfmetrics.registerFont(TTFont('Arial-Bold','/System/Library/Fonts/Supplemental/Arial Bold.ttf'));pdfmetrics.registerFontFamily('Arial',normal='Arial',bold='Arial-Bold',italic='Arial',boldItalic='Arial-Bold')
x=json.loads((W/'content.json').read_text());fusion=json.loads((W/'fusion-status.json').read_text());blocks=x['blocks'][:-1]+[fusion,x['blocks'][-1]]
styles={'title':ParagraphStyle('title',fontName='Arial-Bold',fontSize=17,leading=20,spaceAfter=5,textColor=colors.black),'meta':ParagraphStyle('meta',fontName='Arial',fontSize=8.4,leading=10.5,spaceAfter=10,textColor=colors.black),'body':ParagraphStyle('body',fontName='Arial',fontSize=10.5,leading=13.2,spaceAfter=9,textColor=colors.black),'foot':ParagraphStyle('foot',fontName='Arial',fontSize=8.0,leading=10,spaceBefore=1,textColor=colors.black)}
pdf=O/'Caris-DNA-RNA-research-handoff-TN26-279853-2026-10-05-v2.pdf';story=[Paragraph(x['title'],styles['title']),Paragraph(x['meta'],styles['meta'])]
for b in blocks:story.append(Paragraph('<b>'+b['lead']+'.</b> '+b['body'],styles['body']))
story.append(Paragraph(x['footer'],styles['foot']));doc=SimpleDocTemplate(str(pdf),pagesize=A4,rightMargin=43,leftMargin=43,topMargin=33,bottomMargin=33,title=x['title'],author='Cancer care research',subject='Updated raw Caris DNA/RNA handoff, 5 October 2026');doc.build(story)
r=PdfReader(pdf);assert len(r.pages)==1,f'Expected one page; got {len(r.pages)}'
t='\n\n'.join([x['title'],x['meta']]+[b['lead']+'. '+re.sub('<[^>]+>','',b['body']).replace('&gt;','>') for b in blocks]+[x['footer']]);pdf.with_suffix('.txt').write_text(t+'\n')
for p in [pdf,pdf.with_suffix('.txt'),W/'content.json',W/'fusion-status.json',W/'render.py']:shutil.copy2(p,D/p.name)
for p in W.glob('*questions*.txt'):shutil.copy2(p,D/p.name);shutil.copy2(p,ROOT/'outputs/caris-followup'/p.name)
files=sorted(p for p in D.iterdir() if p.is_file() and p.name!='checksums.sha256');(D/'checksums.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
print(pdf);print('Pages',len(r.pages),'Words',len(t.split()));print('Mirror',D)

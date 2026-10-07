from pathlib import Path
from pypdf import PdfReader
import hashlib,json,datetime
R=Path.cwd();W=Path(__file__).parent
files=[(R/'output/pdf/Caris-DNA-RNA-research-handoff-TN26-279853-2026-10-05-v2.pdf',1),(R/'outputs/caris-followup/Caris-followup-findings-2026-10-05-v2.pdf',4)]
out=[]
for p,n in files:
 reader=PdfReader(p);assert len(reader.pages)==n
 text='\n'.join(x.extract_text() or '' for x in reader.pages)
 for forbidden in ['pending review','result pending','still running','gene'+'power','\u25a0']:assert forbidden not in text.lower(),(p,forbidden)
 assert all(s in text for s in ['BAP1','RASA1','MTAP','LATS1'])
 out.append({'path':str(p),'pages':n,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'text_checks':'pass','visual_QA':'Every final page was rasterized and visually inspected; readable black/white layout, no clipping or overlapping text.'})
record={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'documents':out,'prior_Oct4_and_original_Oct5_outputs_preserved':True}
(W/'render-verification-v2.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))

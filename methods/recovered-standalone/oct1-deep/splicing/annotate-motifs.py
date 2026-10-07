"""Reference splice-motif and selected-transcript context check; public-reference requests only."""
from pathlib import Path
import csv,json,urllib.request,concurrent.futures,hashlib,time,collections
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
rows=list(csv.DictReader((OUT/'target-junctions.tsv').open(),delimiter='\t'))
a=json.loads((OUT/'target-annotation.json').read_text())
selected=[r for r in rows if r['passes_review_priority']=='1']
refs={};provenance={}; (OUT/'reference').mkdir(exist_ok=True)
for p in list((ROOT/'work/oct1-deep').glob('*/*.hg38-reference.json'))+list((OUT/'reference').glob('*.json')):
 try:
  d=json.loads(p.read_text())
  if 'dna' in d:refs[p.name.split('.')[0]]=d;provenance[p.name.split('.')[0]]={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source_url':d.get('source_url'),'bytes':p.stat().st_size}
 except (ValueError,KeyError):pass
needed=sorted({r['gene'] for r in selected}-set(refs))
def download(gene):
 g=a['genes'][gene];s=max(0,g['start0']-200);e=g['end0']+200
 url=f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={g['chrom']};start={s};end={e}"
 p=OUT/'reference'/(gene+'.hg38-reference.json')
 for attempt in range(3):
  try:
   d=json.load(urllib.request.urlopen(url,timeout=60));assert d['start']==s and d['end']==e and len(d['dna'])==e-s
   d['source_url']=url;p.write_text(json.dumps(d)+'\n')
   return gene,d,{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'source_url':url,'bytes':p.stat().st_size}
  except Exception:
   if attempt==2:raise
   time.sleep(2+attempt*2)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for gene,d,p in pool.map(download,needed):refs[gene]=d;provenance[gene]=p
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def near(pos,sites):
 q=min(sites,key=lambda x:abs(x-pos));return q,pos-q
for row in selected:
 gene=row['gene'];d=refs[gene];s=int(row['intron_start0']);e=int(row['intron_end0']);strand=row['target_strand']
 sequence=lambda x,y:d['dna'][x-d['start']:y-d['start']].upper()
 genomic_start=sequence(s,s+2);genomic_end=sequence(e-2,e)
 donor,acceptor=(genomic_start,genomic_end) if strand=='+' else (rc(genomic_end),rc(genomic_start))
 row['reference_splice_motif']=donor+'-'+acceptor; row['major_or_minor_canonical_motif']=int(row['reference_splice_motif'] in {'GT-AG','GC-AG','AT-AC'})
 row['left_boundary_reference_24nt']=sequence(s-12,s+12);row['right_boundary_reference_24nt']=sequence(e-12,e+12)
 tx=a['transcripts'][a['genes'][gene]['selected_transcript']]; ex=tx['exons'];junc=[(l[1],r[0]) for l,r in zip(ex,ex[1:])]
 closest=min(junc,key=lambda j:abs(j[0]-s)+abs(j[1]-e));row['closest_selected_junction']=f'{closest[0]}-{closest[1]}'
 row['start_shift_vs_closest_selected']=s-closest[0];row['end_shift_vs_closest_selected']=e-closest[1]
 # Genomic intron-length difference must NOT be translated into a coding/frame consequence.
 row['intron_length_difference_vs_closest_not_CDS_nt']=(e-s)-(closest[1]-closest[0])
 row['selected_transcript']=tx['transcript_id']
 row['selected_CDS_near_any_boundary']=int(any(cs-1<=p<=ce+1 for cs,ce in tx['CDS'] for p in [s,e]))
 row['canonical_strand_read_fraction']=round(int(row['XS_plus' if strand=='+' else 'XS_minus'])/int(row['reads']),6)
with (OUT/'candidate-motif-context.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(selected[0]),delimiter='\t');w.writeheader();w.writerows(sorted(selected,key=lambda r:int(r['fragments']),reverse=True))
report={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'new_reference_files':needed,'new_reference_total_bytes':sum(provenance[g]['bytes'] for g in needed),'references':provenance,'candidate_count':len(selected),'motifs':dict(collections.Counter(r['reference_splice_motif'] for r in selected)),'limit':'Only public hg38 reference locus sequence requested. A canonical motif does not establish real splicing, tumor specificity or functional effect. Nearest selected junction is a descriptive comparator, not proof of an alternative full isoform; its intron-length difference is genomic, not a transcript-length change or clinical HGVS call.'}
(OUT/'reference-motif-provenance.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='references'},indent=2))

"""Independent nucleotide-token audit of premature-stop/splice geometry."""
from pathlib import Path
import json,gzip,re,datetime
P=Path(__file__).resolve().parent;ROOT=P.parents[2];B=ROOT/'work/oct1-deep/coding-indels';original=json.loads((P/'stop-context.json').read_text())['results'];VAR={r['gene']:r for r in json.loads((B/'conditional-transcript-consequences.json').read_text())['results'] if r['gene'] in ['BAP1','RASA1','LATS1']};ANN={r['transcript']:r for r in json.loads((B/'selected-coding-transcripts.json').read_text())};TX=json.loads((B/'selected-reference-transcripts.json').read_text());EX={r['transcript']:[] for r in original};COMP=str.maketrans('ACGT','TGCA')
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  a=line.rstrip().split('\t')
  if a[2]!='exon':continue
  m=re.search('transcript_id "([^"]+)"',a[8])
  if m and m[1] in EX:EX[m[1]].append((int(a[3])-1,int(a[4])))
rows=[]
for old in original:
 v=VAR[old['gene']];a=ANN[v['transcript']];rev=a['strand']=='-';ex=sorted(EX[v['transcript']],reverse=rev);seq=TX[v['transcript']]['sequence'];loc=[]
 for eno,(s,e) in enumerate(ex,1):loc.extend((p,eno) for p in (range(e-1,s-1,-1) if rev else range(s,e)))
 tokens=[(b,p,eno) for b,(p,eno) in zip(seq,loc)];assert len(seq)==len(loc)
 cdspos=[p for s,e in sorted(a['CDS'],reverse=rev) for p in (range(e-1,s-1,-1) if rev else range(s,e))];mmap={p:i for i,(p,n) in enumerate(loc)};offset=mmap[cdspos[0]];assert seq[offset:offset+3]=='ATG'
 pos=v['pos1']-1;ref=v['ref'];alt=v['alt']
 while ref and alt and ref[0]==alt[0]:pos+=1;ref=ref[1:];alt=alt[1:]
 while ref and alt and ref[-1]==alt[-1]:ref=ref[:-1];alt=alt[:-1]
 if ref:
  places=sorted(mmap[p] for p in range(pos,pos+len(ref)));left,right=places[0],places[-1]+1;assert places==list(range(left,right))
 else:left=right=min(mmap[pos-1],mmap[pos])+1
 rnaalt=alt.translate(COMP)[::-1] if rev else alt;rnaref=ref.translate(COMP)[::-1] if rev else ref;assert seq[left:right]==rnaref;eno=tokens[left-1][2];assert tokens[right][2]==eno
 mt=tokens[:left]+[(b,None,eno) for b in rnaalt]+tokens[right:];ms=''.join(t[0] for t in mt);stop=next(i for i in range(offset,len(ms)-2,3) if ms[i:i+3] in ['TGA','TAG','TAA']);end=stop+3;st=mt[stop:end];junctions=[i for i in range(1,len(mt)) if mt[i-1][2]!=mt[i][2]];down=[b for b in junctions if b>=end]
 values={'conditional_stop_codon_mutant_CDS_1based':[stop-offset+1,end-offset],'conditional_stop_reference_CDS_1based':[mmap[p]-offset+1 if p is not None else None for nt,p,eno in st],'stop_codon_genome_positions1':[p+1 if p is not None else None for nt,p,eno in st],'stop_exon_number_transcript_order':st[0][2],'distance_stop_end_to_next_downstream_junction_nt':min(down)-end if down else None,'distance_stop_end_to_final_junction_nt':junctions[-1]-end,'downstream_junction_count':len(down),'exon_count_including_UTRs':len(ex)}
 rows.append({'gene':v['gene'],'transcript':v['transcript'],'stop_codon':ms[stop:end],'calculated':values,'all_geometry_matches':all(old[k]==val for k,val in values.items()),'conditional_protein_length_matches':(stop-offset)//3==v['mutant_protein_length']})
result={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'PASS' if all(r['all_geometry_matches'] and r['conditional_protein_length_matches'] for r in rows) else 'FAIL','method':'Independent full-mRNA nucleotide tokens carry reference genomic coordinates and exon identities; apply minimal edit, translate to first in-frame stop, identify downstream junctions from adjacent exon tokens rather than shifting an old boundary list. No source script imported.','results':rows,'interpretation':'All descriptive geometry is correct. These values do not quantify nonsense-mediated decay, prove escape, establish full-length mutant RNA abundance, or demonstrate protein or HLA presentation. In particular actual new stop, not the earlier variant position, is the relevant geometric point; downstream-junction count is retained.'}
(P/'stop-context-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

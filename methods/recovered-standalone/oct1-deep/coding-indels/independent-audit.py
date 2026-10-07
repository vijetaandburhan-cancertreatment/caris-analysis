"""Independent bounded audit; no source scripts imported, no network use.
All coordinates are GRCh38, 0-based half-open except VCF positions explicitly pos1.
"""
from pathlib import Path
import json,gzip,re,hashlib,datetime
from collections import defaultdict,Counter
import pysam
P=Path(__file__).resolve().parent; ROOT=P.parents[2]
RAW=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
TARGET={'ANKDD1B','ANKMY1','CEP250','NLRC3','FN3KRP','SLC27A3'}
V=json.loads((P/'expressed-frameshift-shortlist.json').read_text())['candidates']
ANN={r['transcript']:r for r in json.loads((P/'selected-coding-transcripts.json').read_text())}
TX=json.loads((P/'selected-reference-transcripts.json').read_text());PROT=json.loads((P/'selected-reference-proteins.json').read_text())
CONS=json.loads((P/'conditional-transcript-consequences.json').read_text())['results']
EX=defaultdict(list); START=defaultdict(list)
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  a=line.rstrip().split('\t');m=re.search('transcript_id "([^"]+)"',a[8])
  if m and m[1] in TX:
   if a[2]=='exon':EX[m[1]].append((int(a[3])-1,int(a[4])))
   if a[2]=='start_codon':START[m[1]].extend(range(int(a[3])-1,int(a[4])))
COMP=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(COMP)[::-1]
# NCBI standard genetic code, expressed as a fixed table independently of Biopython.
BASES='TCAG';AA='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE=dict(zip((a+b+c for a in BASES for b in BASES for c in BASES),AA))
def trans(s,stop=True):
 a=''.join(CODE.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
 return a.split('*')[0] if stop else a
translation=[]
for prior in CONS:
 tx=prior['transcript'];ann=ANN[tx];var=next(v for v in V if v['chrom']==prior['chrom'] and v['pos1']==prior['pos1'] and v['ref']==prior['ref'] and v['alt']==prior['alt'])
 ex=sorted(EX[tx]);positions=[p for s,e in ex for p in range(s,e)]
 if ann['strand']=='-':positions.reverse()
 genomic_index={p:i for i,p in enumerate(positions)};seq=TX[tx]['sequence'];assert len(seq)==len(positions)
 starts=sorted(genomic_index[p] for p in START[tx]);assert len(starts)==3 and seq[starts[0]:starts[0]+3]=='ATG'
 start=starts[0];refprotein=trans(seq[start:]);row={'gene':prior['gene'],'transcript':tx,'reference_protein_match':refprotein==PROT[tx]['protein'],'prior_status':prior['status']}
 # Edit minimal allele into full reference mRNA (not the prior CDS index).
 pos=var['pos1']-1;ref,alt=var['ref'],var['alt']
 while ref and alt and ref[0]==alt[0]:pos+=1;ref=ref[1:];alt=alt[1:]
 while ref and alt and ref[-1]==alt[-1]:ref=ref[:-1];alt=alt[:-1]
 cds_positions={p for s,e in ann['CDS'] for p in range(s,e)}
 if ref and not all(p in cds_positions for p in range(pos,pos+len(ref))):row['independent_status']='not_wholly_coding'
 elif not ref and (pos not in cds_positions or pos-1 not in cds_positions):row['independent_status']='coding_boundary'
 else:
  if ref:
   loc=sorted(genomic_index[p] for p in range(pos,pos+len(ref)));left,right=loc[0],loc[-1]+1
   if loc!=list(range(left,right)):raise ValueError('Non-contiguous transcript edit')
  else:
   ends=sorted([genomic_index[pos-1],genomic_index[pos]]);assert ends[1]-ends[0]==1;left=right=ends[1]
  rt=ref if ann['strand']=='+' else rc(ref);at=alt if ann['strand']=='+' else rc(alt)
  assert seq[left:right]==rt
  mut=trans((seq[:left]+at+seq[right:])[start:]);pref=0
  while pref<min(len(refprotein),len(mut)) and refprotein[pref]==mut[pref]:pref+=1
  row.update(independent_status='conditional_single_edit',mutant_protein_match=mut==prior.get('mutant_protein'),first_changed_or_lost_residue1=pref+1,mutant_length=len(mut),new_tail_length=len(mut)-pref,net_coding_base_change=len(at)-len(rt))
 translation.append(row)

def window(v):
 p=v['pos1']-1;return json.loads((P/'reference-windows'/f'{v["chrom"]}-{p-80}-{p+len(v["ref"])+80}.json').read_text())
def edit_reference(w,pos,ref,alt):
 i=pos-w['start'];s=w['dna'].upper()
 if i<0 or i+len(ref)>len(s) or s[i:i+len(ref)]!=ref:return None
 return s[:i]+alt+s[i+len(ref):]
# Entire VCF is only ~3k records; no variant caller or BAM rescanning.
vcf=defaultdict(list);vcf_record_count=0
with (RAW/'DNA_TN26-279853.vcf').open() as f:
 for line in f:
  if line.startswith('#'):continue
  vcf_record_count+=1
  a=line.rstrip().split('\t');info=dict(x.split('=',1) if '=' in x else (x,True) for x in a[7].split(';'))
  for alt in a[4].split(','):vcf[a[0]].append({'pos1':int(a[1]),'ref':a[3],'alt':alt,'filter':a[6],'info':info,'format':dict(zip(a[8].split(':'),a[9].split(':')))})
vcf_audit=[]
for v in V:
 w=window(v);target=edit_reference(w,v['pos1']-1,v['ref'],v['alt']);hits=[];near=[]
 for q in vcf[v['chrom']]:
  if abs(q['pos1']-v['pos1'])>70:continue
  altseq=edit_reference(w,q['pos1']-1,q['ref'],q['alt'])
  (hits if altseq==target else near).append(q)
 vcf_audit.append({'gene':';'.join(v['genes']),'pos1':v['pos1'],'matching_alleles':hits,'other_VCF_records_within70bp':near,'match_status':'exact_local_haplotype_present' if hits else 'not_present_as_single_exported_VCF_allele'})
# Independently re-create each public mapping's alternate reference string and
# verify every cited global frequency is linked to exactly the matching allele.
pub=json.loads((P/'expressed-indel-public-context.json').read_text())['results'];public=[]
for r in pub:
 v=next(v for v in V if v['chrom']==r['chrom'] and v['pos1']==r['pos1'] and v['ref']==r['ref'] and v['alt']==r['alt']);w=window(v);s=w['dna'].upper();target=edit_reference(w,v['pos1']-1,v['ref'],v['alt'])
 for hit in r['matches']:
  data=json.loads((P/'public-catalog'/f'{hit["id"]}.variation.json').read_text())['result'];matches=[]
  for m in data.get('mappings',[]):
   if m.get('assembly_name')!='GRCh38' or m.get('seq_region_name')!=v['chrom'][3:] or m.get('strand')!=1:continue
   start,end=m['start']-1,m['end'];i,j=start-w['start'],end-w['start']
   if not 0<=i<=j<=len(s):continue
   for a in m['allele_string'].split('/'):
    alt=a.replace('-','');mut=s[:i]+alt+s[j:]
    if mut==target:matches.append({'start':m['start'],'end':m['end'],'allele':a})
  freq_ok=all(any(z['allele']==f['allele'] for z in matches) and f in data.get('populations',[]) for f in hit['global_frequencies'])
  public.append({'gene':r['gene'],'id':hit['id'],'haplotype_equivalence':bool(matches),'matching_alleles':matches,'cited_frequencies_match_source_and_allele':freq_ok,'maximum_global_frequency':max([float(f['frequency']) for f in hit['global_frequencies']],default=None)})

def extract(read,start,end):
 # Explicit reference-coordinate traversal; both window endpoints must be M/=X.
 rp=read.reference_start;qp=0;sequence=[];quality=[];cover=set();events=[];spans=[]
 for op,n in read.cigartuples:
  if op in (0,7,8):
   s=max(rp,start);e=min(rp+n,end)
   if s<e:
    a=qp+s-rp;b=qp+e-rp;sequence.append(read.query_sequence[a:b]);quality.extend(read.query_qualities[a:b]);cover.update([s,e-1])
   spans.append((rp,rp+n,qp));rp+=n;qp+=n
  elif op==1:
   if start<rp<end:sequence.append(read.query_sequence[qp:qp+n]);quality.extend(read.query_qualities[qp:qp+n])
   events.append(('I',rp,read.query_sequence[qp:qp+n]));qp+=n
  elif op==2:events.append(('D',rp,n));rp+=n
  elif op==3:
   if rp<end and rp+n>start:return None,'splice_inside_window'
   rp+=n
  elif op==4:qp+=n
 if start not in cover or end-1 not in cover:return None,'incomplete'
 if not quality or min(quality)<20:return None,'base_quality'
 return ''.join(sequence),'ok'

def audit_locus(v,bam,flank):
 w=window(v);p=v['pos1']-1;start=p-flank;end=p+len(v['ref'])+flank;s=w['dna'].upper();ref=s[start-w['start']:end-w['start']];alt=ref[:flank]+v['alt']+ref[flank+len(v['ref']):]
 calls=defaultdict(set);info=defaultdict(list);counts=Counter();alt_events=Counter();alt_event_names=defaultdict(set)
 for read in bam.fetch(v['chrom'],start,end):
  if read.flag&(4|256|512|1024|2048) or read.mapping_quality<20:continue
  if read.has_tag('NH') and read.get_tag('NH')>1:counts['NH_multimapper']+=1;continue
  h,reason=extract(read,start,end)
  if h is None:counts[reason]+=1;continue
  call='alt' if h==alt else 'ref' if h==ref else 'other';counts[call]+=1;name=(read.get_tag('RG') if read.has_tag('RG') else '')+'|'+read.query_name;calls[name].add(call)
  if call=='alt' and flank==12:
   row={'name':name,'start0':read.reference_start,'end0':read.reference_end,'reverse':read.is_reverse,'read1':read.is_read1,'mapq':read.mapping_quality,'cigar':read.cigarstring,'NM':read.get_tag('NM') if read.has_tag('NM') else None,'NH':read.get_tag('NH') if read.has_tag('NH') else None,'sequence_window':h,'read_length':read.query_length,'softclip_bases':sum(n for op,n in read.cigartuples if op==4),'mate_start0':read.next_reference_start,'template_length':read.template_length}
   events=[];rp=read.reference_start;qp=0
   for op,n in read.cigartuples:
    if op in (0,7,8):
     for z in range(max(0,w['start']-rp),min(n,w['end']-rp)):
      bp=read.query_sequence[qp+z];refbp=s[rp+z-w['start']]
      if bp!=refbp and read.query_qualities[qp+z]>=20:events.append(('SNV',rp+z,refbp,bp))
     rp+=n;qp+=n
    elif op==1:
     if w['start']<=rp<w['end'] and min(read.query_qualities[qp:qp+n])>=20:events.append(('I',rp,read.query_sequence[qp:qp+n]))
     qp+=n
    elif op==2:
     if w['start']<=rp<w['end']:events.append(('D',rp,n))
     rp+=n
    elif op==3:rp+=n
    elif op==4:qp+=n
   row['local_other_or_target_events']=events;info[name].append(row)
   for event in events:alt_event_names[str(event)].add(name)
 fr=Counter(next(iter(c)) if len(c)==1 else 'discordant' for c in calls.values())
 result={'flank':flank,'reads':dict(counts),'paired_names':dict(fr)}
 if flank==12:
  result['alternate_read_evidence']=[q for name,lst in info.items() if calls[name]=={'alt'} for q in lst]
  result['alternate_pair_endpoint_signatures']=len({tuple(sorted((q['start0'],q['end0'],q['reverse'],q['mate_start0']) for q in lst)) for name,lst in info.items() if calls[name]=={'alt'}})
  result['cooccurring_events_alternate_pair_names']={event:len(names) for event,names in alt_event_names.items()}
 return result
bam_audit=[]
for kind in ['DNA','RNA']:
 with pysam.AlignmentFile(str(RAW/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai')) as bam:
  for v in V:
   if not set(v['genes'])&TARGET:continue
   row={'gene':';'.join(v['genes']),'kind':kind,'chrom':v['chrom'],'pos1':v['pos1'],'ref':v['ref'],'alt':v['alt'],'window_checks':[audit_locus(v,bam,n) for n in (12,20,30,50)]};bam_audit.append(row)
obj={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'VCF_record_count':vcf_record_count,'method':f'Independent fixed-genetic-code full-mRNA edit and stop translation; independent full-window public mapping and frequency equivalence; full {vcf_record_count}-record VCF local haplotype comparisons; bounded six-locus BAM reads using 12/20/30/50-reference-base exact flanks, MAPQ20/BQ20, primary/nonduplicate/QCpass, NH1 when available, paired names collapsed. No source analysis script imported. Counts are not UMI molecules; no matched normal or clinical assertions.','translations':translation,'public_catalog_audit':public,'VCF_audit':vcf_audit,'six_locus_BAM_audit':bam_audit}
(P/'audit-data.json').write_text(json.dumps(obj,indent=2)+'\n')
print('Translation rows',len(translation),'valid',sum(x.get('mutant_protein_match',False) for x in translation),'boundary',sum('mutant_protein_match' not in x for x in translation))
print('Public rows',len(public),'allpass',all(x['haplotype_equivalence'] and x['cited_frequencies_match_source_and_allele'] for x in public))
print('VCF exact equivalence',[(x['gene'],[(h['pos1'],h['filter']) for h in x['matching_alleles']]) for x in vcf_audit])
for r in bam_audit:print(r['gene'],r['kind'],[(x['flank'],x['paired_names']) for x in r['window_checks']])

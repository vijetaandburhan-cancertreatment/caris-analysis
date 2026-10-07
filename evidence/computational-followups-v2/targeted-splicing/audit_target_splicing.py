#!/usr/bin/env python3
"""Bounded, research-only comparison of existing and new four-gene RNA alignments."""
import argparse, collections, csv, hashlib, json, pathlib, re, time
import pysam
GENES={'BAP1','LATS1','LATS2','RASA1'}
EXCLUDE=4|256|512|1024|2048
# Forward-genome edits, zero-based interbase coordinates; insertion has start=end.
EDITS={'BAP1':(52408550,52408561,''),'RASA1':(87332560,87332561,''),'LATS1':(149695191,149695191,'A')}

def attrs(s):return dict(re.findall(r'(\w+) "([^"]*)"',s))

def build_annotation(gtf):
 genes={};tx={}
 with open(gtf) as f:
  for ln in f:
   if ln.startswith('#'):continue
   a=ln.rstrip().split('\t')
   if a[2] not in ('gene','transcript','exon'):continue
   d=attrs(a[8]);g=d.get('gene_name')
   if g not in GENES:continue
   if a[2]=='gene':genes[g]={'chrom':a[0],'start0':int(a[3])-1,'end0':int(a[4]),'strand':a[6]}
   if 'transcript_id' not in d:continue
   t=tx.setdefault(d['transcript_id'],{'gene':g,'exons':[],'tags':[]})
   t['tags'].extend(re.findall(r'tag "([^"]+)"',a[8]))
   if a[2]=='exon':t['exons'].append((int(a[3])-1,int(a[4])))
 known=collections.defaultdict(lambda:collections.defaultdict(set))
 for name,t in tx.items():
  ex=sorted(set(t['exons']))
  for left,right in zip(ex,ex[1:]):known[t['gene']][(left[1],right[0])].add(name)
 return genes,tx,known

def exact_haplotype(r,start,end,ref,alt,flank=10):
 """Query sequence between flanking reference bases, including inserted bases. No soft-clipped reconstruction."""
 l=start-flank;u=end+flank
 mapping={rp:qp for qp,rp in r.get_aligned_pairs(matches_only=False) if rp is not None and qp is not None}
 if l not in mapping or u-1 not in mapping:return None
 q1=mapping[l];q2=mapping[u-1]+1
 if q2<=q1 or r.query_qualities is None or min(r.query_qualities[q1:q2])<20:return None
 seq=r.query_sequence[q1:q2]
 if seq==ref:return 'reference'
 if seq==alt:return 'alternate'
 return 'other'

def qualified_junctions(r,reference,gstart,anchor=12):
 pos=r.reference_start;q=0;out=[];c=r.cigartuples or[];seq=r.query_sequence;qual=r.query_qualities
 for i,(op,n) in enumerate(c):
  if op==3 and i>0 and i+1<len(c) and c[i-1][0] in (0,7,8) and c[i+1][0] in (0,7,8) and min(c[i-1][1],c[i+1][1])>=anchor and qual is not None:
   a=pos-anchor-gstart;b=pos+n-gstart
   if a>=0 and b+anchor<=len(reference) and min(qual[q-anchor:q+anchor])>=20 and seq[q-anchor:q]==reference[a:a+anchor] and seq[q:q+anchor]==reference[b:b+anchor]:out.append((pos,pos+n))
  if op in (0,2,3,7,8):pos+=n
  if op in (0,1,4,7,8):q+=n
 return out

def process(args):
 genes,tx,known=build_annotation(args.gtf);fasta=pysam.FastaFile(args.fasta)
 refseq={g:fasta.fetch(d['chrom'],d['start0'],d['end0']).upper() for g,d in genes.items()}
 contexts={}
 for g,(s,e,ins) in EDITS.items():
  ch=genes[g]['chrom'];ref=fasta.fetch(ch,s-args.flank,e+args.flank).upper();alt=ref[:args.flank]+ins+ref[args.flank+(e-s):]
  contexts[g]={'start0':s,'end0':e,'reference':ref,'alternate':alt}
 frags={g:{} for g in genes};counts=collections.Counter()
 b=pysam.AlignmentFile(args.bam,'rb',index_filename=args.index if args.index else None)
 def records():
  if b.has_index():
   for g,d in genes.items():
    for r in b.fetch(d['chrom'],d['start0'],d['end0']):yield g,r
  else:
   for r in b.fetch(until_eof=True):
    if r.is_unmapped:continue
    for g,d in genes.items():
     if r.reference_name==d['chrom'] and r.reference_start<d['end0'] and r.reference_end>d['start0']:yield g,r
 for g,r in records():
  counts[g+'.all_alignments']+=1
  if r.flag&EXCLUDE or r.mapping_quality<20 or (r.has_tag('NH') and r.get_tag('NH')!=1):continue
  counts[g+'.qualified_alignments']+=1
  key=(r.get_tag('RG') if r.has_tag('RG') else '')+'|'+r.query_name
  f=frags[g].setdefault(key,{'alleles':set(),'junctions':set(),'direct_alt_junctions':set(),'patterns':set()})
  js=qualified_junctions(r,refseq[g],genes[g]['start0']);f['junctions'].update(js)
  f['patterns'].add((r.reference_start,r.cigarstring,r.flag,r.next_reference_start))
  if g in contexts:
   v=contexts[g];allele=exact_haplotype(r,v['start0'],v['end0'],v['reference'],v['alternate'],args.flank)
   if allele:f['alleles'].add(allele)
   if allele=='alternate':f['direct_alt_junctions'].update(js)
 b.close();fasta.close()
 summaries=[];jrows=[];support=[]
 for g in genes:
  ac=collections.Counter();pairs=collections.defaultdict(set);direct=collections.defaultdict(set);alljs=collections.defaultdict(set);patterns=collections.defaultdict(set)
  for key,f in frags[g].items():
   state=next(iter(f['alleles'])) if len(f['alleles'])==1 else ('conflict' if f['alleles'] else 'uncovered')
   ac[state]+=1
   if state in ('alternate','reference','other','conflict'):
    support.append({'gene':g,'fragment':key,'state':state,'junctions':sorted(f['junctions']),'patterns':sorted(f['patterns'])})
   for j in f['junctions']:
    alljs[j].add(key);patterns[j].add(tuple(sorted(f['patterns'])))
    if state=='alternate':pairs[j].add(key)
   if state=='alternate':
    for j in f['direct_alt_junctions']:direct[j].add(key)
  summaries.append({'gene':g,'qualified_query_name_fragments':len(frags[g]),'allele_states':dict(ac),'junction_count':len(alljs),'distinct_start_cigar_patterns_not_UMIs':len({tuple(sorted(f['patterns'])) for f in frags[g].values()})})
  for j,names in sorted(alljs.items()):
   jrows.append({'gene':g,'chrom':genes[g]['chrom'],'intron_start0':j[0],'intron_end0':j[1],'reference_annotation_transcripts':','.join(sorted(known[g].get(j,[]))),'query_name_fragments':len(names),'distinct_fragment_alignment_patterns':len(patterns[j]),'alternate_linked_pairs':len(pairs[j]),'alternate_and_junction_same_read':len(direct[j])})
 out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
 fields=['gene','chrom','intron_start0','intron_end0','reference_annotation_transcripts','query_name_fragments','distinct_fragment_alignment_patterns','alternate_linked_pairs','alternate_and_junction_same_read']
 with (out/'junctions.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(jrows)
 (out/'allele-fragment-evidence.json').write_text(json.dumps(support,indent=2)+'\n')
 report={'bam':str(pathlib.Path(args.bam).resolve()),'label':args.label,'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'pysam_version':pysam.__version__,'filters':{'exclude_flags':EXCLUDE,'MAPQ_min':20,'NH_if_present':1,'base_quality':20,'exact_reference_junction_anchor_each_side':12,'exact_haplotype_flank_each_side':args.flank},'method_notes':['Names group read pairs within RG; not unique original molecules. Distinct alignment patterns are not UMIs.','Duplicate flags excluded, but absence of a duplicate flag does not establish uniqueness; original and new alignment duplicates may differ.','Unstranded bulk tissue RNA cannot identify malignant cell origin.','An annotated junction or short fragment cannot establish a complete expressed isoform, mutant protein production, NMD efficiency or therapeutic sensitivity.','These strict reference-exact flanks can exclude reads with neighboring variants or errors; compare only equal settings.','Insertion/deletion queried by full local haplotype, allowing equivalent CIGAR representations.'],'gene_intervals':genes,'variant_contexts':contexts,'alignment_counts':dict(counts),'summaries':summaries}
 (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'label':args.label,'summaries':summaries},indent=2))

def main():
 p=argparse.ArgumentParser();p.add_argument('--bam',required=True);p.add_argument('--index');p.add_argument('--gtf',required=True);p.add_argument('--fasta',required=True);p.add_argument('--out',required=True);p.add_argument('--label',required=True);p.add_argument('--flank',type=int,default=10);a=p.parse_args();process(a)
if __name__=='__main__':main()

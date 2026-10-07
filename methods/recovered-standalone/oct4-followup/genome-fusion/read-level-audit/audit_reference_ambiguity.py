#!/usr/bin/env python3
"""Stream local reference sequences for exact junction/anchor alternatives."""
import argparse,collections,gzip,hashlib,json,pathlib,resource,time
import ahocorasick_rs,pysam

def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]

def scan_genome(path,auto,patterns,maxlen,chunk_bases=262144,examples=30):
 counts=collections.Counter();hits=collections.defaultdict(list);basecount=0;records=0
 def run_chunk(contig,chunk,tail,consumed):
  joined=tail+chunk;origin=consumed-len(tail)
  for i,a,b in auto.find_matches_as_indexes(joined,overlapping=True):
   if b<=len(tail):continue
   counts[i]+=1
   if len(hits[i])<examples:hits[i].append({'contig':contig,'start0':origin+a,'end0':origin+b})
  return joined[-maxlen+1:] if maxlen>1 else '',consumed+len(chunk)
 with open(path) as f:
  contig=None;tail='';consumed=0;buf=[];nb=0
  for line in f:
   if line.startswith('>'):
    if buf:tail,consumed=run_chunk(contig,''.join(buf),tail,consumed);buf=[];nb=0
    contig=line[1:].split()[0];records+=1;tail='';consumed=0
   else:
    s=line.strip().upper();buf.append(s);nb+=len(s);basecount+=len(s)
    if nb>=chunk_bases:tail,consumed=run_chunk(contig,''.join(buf),tail,consumed);buf=[];nb=0
  if buf:run_chunk(contig,''.join(buf),tail,consumed)
 return {'total_reference_bases':basecount,'contigs':records,'counts':dict(counts),'first_loci':dict(hits),'maximum_saved_loci_per_pattern':examples}

def main():
 p=argparse.ArgumentParser();p.add_argument('--candidates',required=True);p.add_argument('--genome',required=True);p.add_argument('--transcripts',required=True);p.add_argument('--out',required=True);args=p.parse_args();out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True);t0=time.time();rows=json.load(open(args.candidates));index=collections.defaultdict(list);idfreq=collections.Counter(r['candidate_id'] for r in rows)
 for n,r in enumerate(rows):
  rid=r['candidate_id'] if idfreq[r['candidate_id']]==1 else f"{r['candidate_id']}__{r['source_set']}__{n}"
  for ptn in r['junction_patterns']:
   arm=ptn['arm_bases'];s=ptn['sequence']
   for role,word in [('junction',s),('left_anchor',s[:arm]),('right_anchor',s[-arm:])]:
    for seq,orientation in [(word,'forward'),(rc(word),'reverse')]:index[seq].append({'candidate':rid,'kind':role,'arm':arm,'orientation':orientation})
 patterns=sorted(index);auto=ahocorasick_rs.AhoCorasick(patterns) if patterns else None
 if not patterns:(out/'reference-ambiguity.json').write_text(json.dumps({'status':'unassessable','reason':'no exact junction patterns'}));return
 maxlen=max(map(len,patterns))
 # Exact-match/overlap and chunk-coordinate public regression on a synthetic reference, independent of patient sequence.
 test_patterns=['ACGTAC','GTACGT','TTTTG','CGT'];test_auto=ahocorasick_rs.AhoCorasick(test_patterns);testseq='N'*13+'ACGTACGTACGT'+'N'*3+'TTTTG'+'N'*2;test=out/'public-chunk-regression.fa';test.write_text('>public1\n'+'\n'.join(testseq[i:i+5] for i in range(0,len(testseq),5))+'\n>public2\nACGTAC\n')
 small=scan_genome(test,test_auto,test_patterns,max(map(len,test_patterns)),chunk_bases=7)
 expected=collections.Counter();locs=collections.defaultdict(list)
 for contig,s in [('public1',testseq),('public2','ACGTAC')]:
  for i,w in enumerate(test_patterns):
   for a in range(len(s)-len(w)+1):
    if s[a:a+len(w)]==w:expected[i]+=1;locs[i].append({'contig':contig,'start0':a,'end0':a+len(w)})
 assert dict(expected)==small['counts'] and dict(locs)==small['first_loci'];(out/'reference-scan-public-control.json').write_text(json.dumps({'passed':True,'counts':small,'expected_counts':dict(expected)},indent=2))
 genome=scan_genome(args.genome,auto,patterns,maxlen);(out/'genome-pattern-hits.json').write_text(json.dumps(genome,indent=2));print(json.dumps({'stage':'genome_done','bases':genome['total_reference_bases'],'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}),flush=True)
 counts=collections.Counter();records=collections.Counter();genes=collections.defaultdict(collections.Counter);hits=collections.defaultdict(list);total_bases=total_transcripts=0
 with pysam.FastxFile(args.transcripts) as f:
  for r in f:
   parts=r.name.split('|');gene=parts[5] if len(parts)>5 else r.name;transcript=parts[0];s=r.sequence.upper();total_bases+=len(s);total_transcripts+=1;seen=set()
   for i,a,b in auto.find_matches_as_indexes(s,overlapping=True):
    counts[i]+=1;seen.add(i)
    if len(hits[i])<30:hits[i].append({'transcript':transcript,'gene':gene,'start0':a,'end0':b})
   for i in seen:records[i]+=1;genes[i][gene]+=1
 trans={'total_bases':total_bases,'total_transcripts':total_transcripts,'occurrences':dict(counts),'matching_transcripts':dict(records),'matching_genes':{i:dict(v) for i,v in genes.items()},'first_hits':dict(hits),'maximum_saved_examples_per_pattern':30};(out/'transcript-pattern-hits.json').write_text(json.dumps(trans,indent=2))
 result={'status':'complete','genome_path':args.genome,'transcripts_path':args.transcripts,'patterns':[{'sequence':s,'roles':index[s]} for s in patterns],'genome':genome,'transcripts':trans,'public_chunk_coordinate_and_overlap_control':True,'candidate_source_sha256':hashlib.sha256(pathlib.Path(args.candidates).read_bytes()).hexdigest(),'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'limits':['Exact reference matching is not comprehensive alignment/mappability analysis; mismatched or gapped competing alignments, absent haplotypes and unannotated isoforms remain possible.','An individual anchor may cross a normal exon junction and have no contiguous-genome match; absence is not automatically uniqueness.','Normal-reference transcript junction matches may explain an RNA junction without a tumor-specific DNA rearrangement; interpretation must use longer context/coordinates.','Genomic hit positions are0-based half-open; strand follows the pattern orientation relative to the reported fusion sequence.']};(out/'reference-ambiguity.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k] for k in ['status','elapsed_seconds','peak_RSS_bytes']}),flush=True)
if __name__=='__main__':main()

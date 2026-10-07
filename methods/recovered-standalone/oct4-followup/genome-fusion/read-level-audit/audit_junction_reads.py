#!/usr/bin/env python3
"""Exact raw-read junction audit; no alignment and no actionability inference."""
import argparse,collections,gzip,hashlib,itertools,json,pathlib,resource,time
import ahocorasick_rs,pysam

def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def name(s):return s[:-2] if s.endswith(('/1','/2')) else s

def main():
 p=argparse.ArgumentParser();p.add_argument('--candidates',required=True);p.add_argument('--r1',required=True);p.add_argument('--r2',required=True);p.add_argument('--expected-pairs',type=int,required=True);p.add_argument('--out',required=True);p.add_argument('--chimeric');args=p.parse_args();out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True);t0=time.time()
 rows=json.load(open(args.candidates));idfreq=collections.Counter(r['candidate_id'] for r in rows);sources=collections.defaultdict(list);metadata=[];name_to_candidates=collections.defaultdict(set)
 for i,r in enumerate(rows):
  rid=r['candidate_id'] if idfreq[r['candidate_id']]==1 else f"{r['candidate_id']}__{r['source_set']}__{i}"
  r['audit_key']=rid;metadata.append(r)
  for n in r['read_names']:name_to_candidates[name(n)].add(rid)
  for pat in r['junction_patterns']:
   s=pat['sequence'];arm=pat['arm_bases']
   for seq,ori in [(s,'forward'),(rc(s),'reverse')]:sources[seq].append((rid,arm,ori))
 pats=sorted(sources);auto=ahocorasick_rs.AhoCorasick(pats) if pats else None
 # Algorithm regression before private input. Overlaps and repeated patterns must agree with plain string slicing.
 for seq in [s+s[:10] for s in pats[:100]]+[rc(s) for s in pats[:100]]:
  got=sorted(auto.find_matches_as_indexes(seq,overlapping=True));expected=sorted((i,k,k+len(w)) for i,w in enumerate(pats) for k in range(len(seq)-len(w)+1) if seq[k:k+len(w)]==w);assert got==expected
 stats={r['audit_key']:{'caller_named_pairs_recovered':0,'caller_names_missing':set(r['read_names']),'raw_exact':{str(a):{f'Q{q}':{'pairs':0,'caller_named_pairs':0,'sequence_families':set(),'orientation_normalized_sequence_families':set(),'caller_named_sequence_families':set(),'R1':0,'R2':0,'forward':0,'reverse':0,'both_mates':0,'minimum_marker_end_distance_hist':collections.Counter(),'junction_cycle_hist':collections.Counter()} for q in [20,30]} for a in [15,20,25]}} for r in rows}
 paired_out=[gzip.open(out/f'audit-selected.R{x}.fastq.gz','wt',compresslevel=5) for x in [1,2]];ev=gzip.open(out/'read-evidence.jsonl.gz','wt',compresslevel=5);counts=collections.Counter();selectednames=set();pairnamehashes=set();exception=None
 try:
  with pysam.FastxFile(args.r1) as f1,pysam.FastxFile(args.r2) as f2:
   for r1,r2 in itertools.zip_longest(f1,f2):
    if r1 is None or r2 is None:raise ValueError('unequal FASTQ EOF')
    if name(r1.name)!=name(r2.name):raise ValueError('mate-name disagreement')
    counts['pairs_scanned']+=1
    if counts['pairs_scanned']%1000000==0:
     pr={'counts':dict(counts),'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(out/'progress.json').write_text(json.dumps(pr,indent=2));print(json.dumps(pr),flush=True)
     if pr['peak_RSS_bytes']>2500000000:raise RuntimeError('RSS guard exceeded2.5GB')
    n=name(r1.name);reported=name_to_candidates.get(n,set());hits=[]
    for mate,r in enumerate([r1,r2]):
     if r.quality is None or len(r.sequence)!=len(r.quality):raise ValueError('sequence/quality mismatch')
     if auto:
      for idx,start,end in auto.find_matches_as_indexes(r.sequence,overlapping=True):
       q=min(ord(c)-33 for c in r.quality[start:end]);edge=min(start,len(r.sequence)-end)
       for rid,arm,orientation in sources[pats[idx]]:hits.append({'candidate':rid,'arm':arm,'mate':mate+1,'start0':start,'end0':end,'min_base_quality':q,'minimum_end_distance':edge,'junction_cycle1':start+arm+1,'orientation':orientation})
    if not hits and not reported:continue
    counts['selected_pairs']+=1;selectednames.add(n);nh=hashlib.sha256(n.encode()).hexdigest()
    if nh in pairnamehashes:raise ValueError('duplicate selected pair name')
    pairnamehashes.add(nh)
    for h,r in zip(paired_out,[r1,r2]):h.write('@'+r.name+(' '+r.comment if r.comment else '')+'\n'+r.sequence+'\n+\n'+r.quality+'\n')
    ev.write(json.dumps({'name_sha256':nh,'reported_candidates':sorted(reported),'hits':hits})+'\n')
    seqpair=r1.sequence+'|'+r2.sequence;seqhash=hashlib.sha256(seqpair.encode()).hexdigest();familyhash=hashlib.sha256(min(seqpair,rc(r2.sequence)+'|'+rc(r1.sequence)).encode()).hexdigest()
    for rid in reported:stats[rid]['caller_named_pairs_recovered']+=1
    for rid in {h['candidate'] for h in hits}:
     for arm in [15,20,25]:
      for q in [20,30]:
       hh=[h for h in hits if h['candidate']==rid and h['arm']==arm and h['min_base_quality']>=q]
       if not hh:continue
       s=stats[rid]['raw_exact'][str(arm)][f'Q{q}'];s['pairs']+=1;s['caller_named_pairs']+=rid in reported;s['sequence_families'].add(seqhash);s['orientation_normalized_sequence_families'].add(familyhash)
       if rid in reported:s['caller_named_sequence_families'].add(seqhash)
       for mate in [1,2]:s[f'R{mate}']+=any(h['mate']==mate for h in hh)
       for ori in ['forward','reverse']:s[ori]+=any(h['orientation']==ori for h in hh)
       s['both_mates']+=len({h['mate'] for h in hh})==2
       # One most-interior hit per pair, deterministic tie-break, for cycle distribution.
       best=max(hh,key=lambda h:(h['minimum_end_distance'],-h['mate'],-h['start0']));s['minimum_marker_end_distance_hist'][str(best['minimum_end_distance'])]+=1;s['junction_cycle_hist'][str(best['junction_cycle1'])]+=1
    if counts['selected_pairs']>200000:raise RuntimeError('More than200,000 selected pairs; investigate common junction patterns before unbounded extraction')
  if counts['pairs_scanned']!=args.expected_pairs:raise ValueError('unexpected total pairs')
 except Exception as exc:exception=repr(exc);raise
 finally:
  for h in paired_out:h.close()
  ev.close()
  if exception:(out/'failure.json').write_text(json.dumps({'error':exception,'counts':dict(counts)}))
 for r in rows:
  rid=r['audit_key'];stats[rid]['caller_names_missing']=sorted(n for n in r['read_names'] if name(n) not in selectednames)
  for aa in stats[rid]['raw_exact'].values():
   for ss in aa.values():
    for k,v in list(ss.items()):
     if isinstance(v,set):ss[k]=len(v)
 junction_count=collections.Counter()
 if args.chimeric:
  with open(args.chimeric) as src,gzip.open(out/'selected-Chimeric.out.junction.gz','wt') as dest:
   for line in src:
    if line.startswith('#'):dest.write(line);continue
    col=line.rstrip('\n').split('\t')
    if len(col)>9 and name(col[9]) in selectednames:
     dest.write(line);junction_count['selected_lines']+=1
     for rid in name_to_candidates.get(name(col[9]),[]):junction_count[rid]+=1
 result={'status':'complete','counts':dict(counts),'input_paths':[args.r1,args.r2],'input_bytes':[pathlib.Path(x).stat().st_size for x in [args.r1,args.r2]],'candidate_source_sha256':hashlib.sha256(pathlib.Path(args.candidates).read_bytes()).hexdigest(),'rows':rows,'statistics':stats,'chimeric_line_counts':dict(junction_count),'public_algorithm_regression':True,'normal_EOF_pair_names_and_quality_lengths':True,'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'method':'Full raw FASTQ exact balanced junction markers, every marker base at Q20/Q30; overlaps enabled; both orientations; count once per query-name pair/arm/quality. Caller-named and all-raw counts separated. Read names and distinct sequence families are not UMI-defined independent molecules.','limits':['Caller-assembled sequence is used to construct markers, so this is raw-read corroboration, not independent biological validation.','No exact marker is not proof of no fusion: sequence mismatch, microhomology, clipped/short anchor or variable splice isoform can prevent exact matches.','Any unlabeled exact match may have an alternate normal/paralog alignment; reference-ambiguity audit must follow.','No tumor-specificity, DNA rearrangement or drug actionability inference from a junction or gene name alone.']}
 (out/'junction-read-audit.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k] for k in ['status','counts','elapsed_seconds','peak_RSS_bytes']}),flush=True)
if __name__=='__main__':main()

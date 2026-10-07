from pathlib import Path
import sys,json,time,resource,gzip,re,hashlib
from collections import defaultdict,Counter
from Bio import SeqIO
ROOT=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-hla-allele-support');sys.path.insert(0,str(ROOT/'deps'));import ahocorasick
REF=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools/arcasHLA/dat');t0=time.time();KS=[31,51,71,91,111,131]
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def two(a):
 m=re.match(r'([A-Z0-9]+\*\d+:\d+)',a);return m.group(1) if m else a
old=json.load(open(ROOT/'markers.json'));targets={t:min(rs,key=lambda x:int(x['id'])) for t,rs in old['targets'].items()};cands={}
for t,r in targets.items():
 s=r['sequence'];other=next(o['sequence'] for tt,o in targets.items() if tt!=t)
 for k in KS:
  for p in range(len(s)-k+1):
   w=s[p:p+k]
   if w in other or rc(w) in other:continue
   if w not in cands:cands[w]={'allele':t,'k':k,'positions':[],'groups':set(),'nontarget_examples':set()}
   cands[w]['positions'].append(p)
keys=list(cands);a=ahocorasick.Automaton();idx=defaultdict(set)
for i,w in enumerate(keys):
 for ww in {w,rc(w)}:idx[ww].add(i)
for w,ii in idx.items():a.add_word(w,tuple(ii))
a.make_automaton();counts=Counter()
for r in SeqIO.parse(REF/'IMGTHLA/hla.dat','imgt'):
 for f in r.features:
  if f.type!='CDS':continue
  al=f.qualifiers.get('allele',[''])[0].removeprefix('HLA-');s=str(f.extract(r.seq)).upper();counts['CDS_records']+=1
  ids={i for p,ii in a.iter(s) for i in ii}
  for i in ids:
   d=cands[keys[i]];d['groups'].add(two(al))
   if two(al)!=d['allele'] and len(d['nontarget_examples'])<5:d['nontarget_examples'].add(al)
for d in cands.values():d['non_HLA_B_transcript_hits']=set()
ref=ROOT/'gencode.v37.transcripts.fa.gz'
with gzip.open(ref,'rt') as f:
 for r in SeqIO.parse(f,'fasta'):
  gene=r.id.split('|')[5];counts['transcripts']+=1
  if gene=='HLA-B':continue
  for i in {i for p,ii in a.iter(str(r.seq).upper()) for i in ii}:cands[keys[i]]['non_HLA_B_transcript_hits'].add(gene)
markers=[]
for w,d in cands.items():
 markers.append({'sequence':w,'allele':d['allele'],'k':d['k'],'positions':[[targets[d['allele']]['id'],p] for p in d['positions']],'compatible_two_field':sorted(d['groups']),'nontarget_examples':sorted(d['nontarget_examples']),'all_reference_exclusive':d['groups']=={d['allele']},'B_gene_exclusive':all(x.split('*')[0]=='B' for x in d['groups']),'passes_transcript_paralog_mask':not d['non_HLA_B_transcript_hits'],'non_HLA_B_transcript_hits':sorted(d['non_HLA_B_transcript_hits']),'pair_discriminating':True})
summary={'status':'prepared','reference':'IMGT3.46 all CDS plus GENCODE37 transcripts; common representative coding sequences B*50:01:01/B*52:01:01','counts':dict(counts),'by_allele_length':{},'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
for t in targets:
 summary['by_allele_length'][t]={str(k):{'pair_discriminating':sum(m['allele']==t and m['k']==k for m in markers),'paralog_safe_pair_discriminating':sum(m['allele']==t and m['k']==k and m['B_gene_exclusive'] and m['passes_transcript_paralog_mask'] for m in markers),'all_reference_exclusive':sum(m['allele']==t and m['k']==k and m['all_reference_exclusive'] and m['passes_transcript_paralog_mask'] for m in markers)} for k in KS}
with gzip.open(ROOT/'common-haplotype-markers.json.gz','wt') as f:json.dump({'summary':summary,'targets':targets,'markers':markers},f)
(ROOT/'common-haplotype-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)

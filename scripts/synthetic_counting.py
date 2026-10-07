"""Historical query-name-collapse parser, isolated for synthetic regression.

The incoming pileup must already have the documented upstream read/BQ/MAPQ
filters applied and the expected nine-column encoding. This parser does not
perform those filters, align reads or validate patient results. Tests use only
invented names/bases; it is not a UMI molecule counter.

Copied function (unchanged) from the authored population-panel counter.
Original source SHA256: a22184bc3b70c9a10ed368963d097e8452486d164c5354c3ee9deac5739fa6cc
Function-source SHA256: bad32721ed47d4274c722e50ce147ad9ea9613942e271cdf5512bc435557c65a
"""
import collections

def parse_line(ln,marker):
 a=ln.rstrip('\n').split('\t')
 # With no FASTA, aligned bases are explicit. Output order from samtools is
 # bases/BQ, MQ, QNAME, BP5 (qualified below using synthetic control).
 if len(a)!=9:raise ValueError(f'Expected9fields,got{len(a)}')
 if int(a[3])==0:
  bases='';quals='';mq='';names=[];bp=[]
 else:
  bases=a[4];quals=a[5];mq=a[6];names=a[7].split(',');bp=[int(x) for x in a[8].split(',')]
 if any(len(x)!=len(bases) for x in [quals,mq,names,bp]):raise ValueError('Pileupparallelfieldlengthmismatch')
 groups=collections.defaultdict(list)
 for b,q,m,n,p in zip(bases,quals,mq,names,bp):
  b=marker['ref'] if b in '.,' else b.upper()
  groups[n].append((b,ord(q)-33,ord(m)-33,p))
 ctr=collections.Counter();read=collections.Counter();edge=collections.Counter();strand=collections.Counter()
 for b in bases:
  if b in 'ACGTacgt':strand[(b.upper(),'rev' if b.islower() else 'fwd')]+=1
 for n,items in groups.items():
  alleles={b for b,q,m,p in items if b in 'ACGT'}
  for b,q,m,p in items:read[b]+=1
  if len(alleles)>1:ctr['discordant']+=1;continue
  if not alleles:ctr['nonACGT']+=1;continue
  b=next(iter(alleles));lab='ref' if b==marker['ref'] else 'alt' if b==marker['alt'] else 'other';ctr[lab]+=1
  if any(p>5 for bb,q,m,p in items if bb==b):edge[lab]+=1
 depth=ctr['ref']+ctr['alt'];af=ctr['alt']/depth if depth else None
 return {**marker,'read_ref':read[marker['ref']],'read_alt':read[marker['alt']],'fragment_ref':ctr['ref'],'fragment_alt':ctr['alt'],'fragment_other':ctr['other'],'fragment_discordant':ctr['discordant'],'fragment_nonACGT':ctr['nonACGT'],'fragment_depth_refalt':depth,'fragment_ALT_fraction':af,'minor_allele_fraction':min(af,1-af) if af is not None else None,'alt_forward_reads':strand[(marker['alt'],'fwd')],'alt_reverse_reads':strand[(marker['alt'],'rev')],'ref_forward_reads':strand[(marker['ref'],'fwd')],'ref_reverse_reads':strand[(marker['ref'],'rev')],'fragment_ref_BP5_gt5':edge['ref'],'fragment_alt_BP5_gt5':edge['alt'],'reported_pileup_depth':int(a[3])}

"""Physical paired-fragment phasing between a public BAP1 SNP and frameshift."""
import pathlib,json,collections,hashlib,statistics
import pysam,py2bit
OUT=pathlib.Path(__file__).resolve().parent;ROOT=OUT.parents[2];SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
CH='chr3';SNP=52408337;DS=52408550;DE=DS+11
tb=py2bit.open(str(OUT/'public/hg38.2bit'));assert tb.sequence(CH,SNP,SNP+1).upper()=='G'
def layout(r):
 p=r.reference_start;q=0;mp={};skips=[]
 for op,n in r.cigartuples:
  if op in [0,7,8]:mp.update((p+i,q+i) for i in range(n));p+=n;q+=n
  elif op in [1,4]:q+=n
  elif op==2:p+=n
  elif op==3:skips.append((p,p+n));p+=n
 return mp,skips
bam=pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai'))
allreads=list(bam.fetch(CH,SNP-250,DE+250));results=[]
for flank in [3,10,20]:
 left=DS-flank;right=DE+flank;wt=tb.sequence(CH,left,right).upper();mt=wt[:flank]+wt[flank+11:];fr=collections.defaultdict(lambda:{'SNP':set(),'deletion':set(),'SNP_clean':set(),'deletion_clean':set(),'alignments':set(),'read1_calls':[],'read2_calls':[]});marginal=collections.Counter()
 for r in allreads:
  if r.flag&(4|8|256|512|1024|2048) or not r.is_proper_pair or r.mapping_quality<30 or r.query_qualities is None:continue
  mp,skips=layout(r);key=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name);d=fr[key];calls={};soft=any(op==4 for op,n in r.cigartuples)
  if SNP in mp:
   q=mp[SNP]
   if r.query_qualities[q]>=25 and r.query_sequence[q] in 'ACGT':
    a=r.query_sequence[q];d['SNP'].add(a);calls['SNP']=a;marginal['SNP_'+a+'_reads']+=1
    if not soft and min(q-r.query_alignment_start,r.query_alignment_end-1-q)>=5:d['SNP_clean'].add(a)
  if left in mp and right-1 in mp and not any(s<right and e>left for s,e in skips):
   qa=mp[left];qb=mp[right-1]+1
   if qb>qa and min(r.query_qualities[qa:qb])>=25:
    obs=r.query_sequence[qa:qb];a='del11' if obs==mt else 'reference' if obs==wt else 'other';d['deletion'].add(a);calls['deletion']=a;marginal['deletion_'+a+'_reads']+=1
    if not soft and min(qa-r.query_alignment_start,r.query_alignment_end-qb)>=5:d['deletion_clean'].add(a)
  if calls:
   d['alignments'].add((r.reference_start,r.reference_end,r.cigarstring,r.is_read1,r.is_reverse,abs(r.template_length)))
   d['read1_calls' if r.is_read1 else 'read2_calls'].append(calls)
 counts=collections.Counter();clean=collections.Counter();details=[];marginals=collections.Counter()
 for key,d in fr.items():
  for locus in ['SNP','deletion']:
   if len(d[locus])==1:marginals[locus+'_'+next(iter(d[locus]))+'_fragments']+=1
   elif len(d[locus])>1:marginals[locus+'_discordant_fragments']+=1
  if not d['SNP'] or not d['deletion']:continue
  if len(d['SNP'])!=1 or len(d['deletion'])!=1:counts['discordant']+=1;continue
  a=next(iter(d['SNP']));b=next(iter(d['deletion']));label=a+'__'+b;counts[label]+=1
  isclean=d['SNP_clean']=={a} and d['deletion_clean']=={b}
  if isclean:clean[label]+=1
  sameread=any('SNP' in x and 'deletion' in x for x in d['read1_calls']+d['read2_calls'])
  details.append({'fragment_identifier_sha256':hashlib.sha256(('|'.join(key)).encode()).hexdigest(),'SNP':a,'deletion':b,'both_loci_clean':isclean,'same_read_both_loci':sameread,'alignment_patterns':[list(z) for z in sorted(d['alignments'])]})
 results.append({'flank_bases':flank,'SNP_pos1':SNP+1,'SNP_ref':'G','SNP_alt':'C','deletion_start0':DS,'deletion_end0':DE,'joint_fragment_counts':dict(counts),'joint_clean_fragment_counts':dict(clean),'marginal_fragment_counts':dict(marginals),'marginal_read_counts':dict(marginal),'joint_fragment_details':details})
(OUT/'BAP1-SNP-frameshift-fragment-phase.json').write_text(json.dumps({'scope':'same-read/same-query-name paired fragment phase, not chromosome-scale/full-transcript phasing','filters':'DNA properpair MAPQ>=30 BQ>=25; exclude unmapped/mateunmapped,secondary,supplementary,QCfail,duplicateflags. Exact localdeletion/reference haplotypes and publicreferenceSNP. Discordantmates excluded. Clean requires no softclips and both locus observations>=5ntfromalignedreadends.','limitations':['Query-name fragments are not UMI unique molecules.','Physical linkage does not by itself establish somaticstatus,loss of oppositeallele,proteinfunction,orclinical biallelicinactivation.','Only fragments coveringbothloci informphase; shortinsert size limitsrepresentativeness.'],'results':results},indent=2)+'\n')
for r in results:print(r['flank_bases'],r['joint_fragment_counts'],r['joint_clean_fragment_counts'])

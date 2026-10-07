"""Bounded candidate generation in22 additional targeted-therapy/DNA-repair genes."""
import pathlib,json,csv,urllib.request,time,collections,statistics
import pysam
from Bio.Seq import Seq
B=pathlib.Path(__file__).resolve().parent;ROOT=B.parents[2];OLD=ROOT/'work/oct1-analysis'
SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
A=json.loads((B/'target-annotation.json').read_text())['genes'];start=time.time()
GENES=['PIK3CA', 'FGFR1', 'FGFR2', 'FGFR3', 'ALK', 'ROS1', 'RET', 'NTRK1', 'NTRK2', 'NTRK3', 'IDH1', 'IDH2', 'KIT', 'PDGFRA', 'MAP2K1', 'MAP2K2', 'BRCA1', 'BRCA2', 'PALB2', 'ATM', 'POLE', 'POLD1']
def merge(a):
 out=[]
 for s,e in sorted(a):
  if out and s<=out[-1][1]:out[-1][1]=max(e,out[-1][1])
  else:out.append([s,e])
 return out
def eligible(r):return not (r.flag & (4|256|512|1024|2048)) and r.mapping_quality>=20
def save(name,x):(B/name).write_text(json.dumps(x,indent=2)+'\n')
refs={}
for gene in GENES:
 g=A[gene];p=B/f'{gene}.hg38-reference.json';url=f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={g['chrom']};start={g['gene_start0']};end={g['gene_end0']}"
 if not p.exists():
  for i in range(3):
   try:
    with urllib.request.urlopen(url,timeout=60) as r:x=json.load(r)
    assert len(x['dna'])==g['gene_end0']-g['gene_start0'];x['source_url']=url;save(p.name,x);break
   except Exception:
    if i==2:raise
    time.sleep(2)
 refs[gene]=json.loads(p.read_text());print('reference',gene,len(refs[gene]['dna']),flush=True)

def annotate(gene,p,ref,alt):
 g=A[gene];seq=refs[gene]['dna'].upper();rs=g['gene_start0'];plus=''.join(seq[s-rs:e-rs] for s,e in g['selected_CDS']);off=0;inside=False
 for s,e in g['selected_CDS']:
  if s<=p<e:off+=p-s;inside=(p+len(ref)<=e);break
  off+=e-s
 if not inside:return {'consequence':'CDS_flank_or_exon_boundary','protein':None,'coding':None}
 assert plus[off:off+len(ref)]==ref
 mut=plus[:off]+alt+plus[off+len(ref):];cds=plus if g['strand']=='+' else str(Seq(plus).reverse_complement());mut=mut if g['strand']=='+' else str(Seq(mut).reverse_complement());rp=str(Seq(cds).translate());mp=str(Seq(mut[:len(mut)//3*3]).translate());assert rp.startswith('M') and '*' not in rp
 diff=next((i for i,(r,a) in enumerate(zip(rp,mp)) if r!=a),None)
 if len(ref)==len(alt)==1:
  c=off if g['strand']=='+' else len(cds)-off-1;coding=f'c.{c+1}{cds[c]}>{mut[c]}'
 else:coding='genomic_indel_not_HGVS_normalized'
 if diff is None:consequence='synonymous';protein='='
 elif mp[diff]=='*':consequence='stop_gain';protein=rp[diff]+str(diff+1)+'*'
 elif len(cds)!=len(mut):consequence='frameshift' if (len(cds)-len(mut))%3 else 'inframe_indel';protein=(rp[diff]+str(diff+1)+mp[diff]+'fs') if consequence=='frameshift' else 'in-frame indel; precise HGVS not assigned by screen'
 else:consequence='missense';protein=rp[diff]+str(diff+1)+mp[diff]
 return {'consequence':consequence,'protein':protein,'coding':coding,'transcript':g['selected_transcript']}

cands=[];summ=[]
with pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(OLD/'DNA_TN26-279853.bam.bai')) as bam:
 for gene in GENES:
  g=A[gene];seq=refs[gene]['dna'].upper();rs=g['gene_start0'];regions=merge([(s-4,e+4) for s,e in g['selected_CDS']]);depth={};indels=collections.defaultdict(set);indread=collections.Counter();indstrands=collections.defaultdict(set)
  for s,e in regions:
   cov=bam.count_coverage(g['chrom'],s,e,quality_threshold=20,read_callback=eligible)
   for i,cs in enumerate(zip(*cov)):
    p=s+i;total=sum(cs);depth[p]=total;rb=seq[p-rs]
    for alt,n in zip('ACGT',cs):
     if alt!=rb and n>=10 and total>=20 and n/total>=.05:
      cands.append({'gene':gene,'chrom':g['chrom'],'pos1':p+1,'ref':rb,'alt':alt,'screen_type':'SNV','screen_alt_reads':n,'screen_anchor_read_depth':total,'screen_fraction':n/total,**annotate(gene,p,rb,alt)})
   for r in bam.fetch(g['chrom'],s,e):
    if not eligible(r) or r.query_qualities is None:continue
    pos=r.reference_start;q=0
    for op,n in r.cigartuples or []:
     if op in (0,7,8):pos+=n;q+=n
     elif op==1:
      if s<=pos-1<e and q>0 and min(r.query_qualities[q-1:q+n])>=20:
       rb=seq[pos-1-rs];alt=rb+r.query_sequence[q:q+n];key=(pos-1,rb,alt);indels[key].add((r.get_tag('RG') if r.has_tag('RG') else '',r.query_name));indread[key]+=1;indstrands[key].add(r.is_reverse)
      q+=n
     elif op==2:
      if s<=pos-1<e and q>0 and r.query_qualities[q-1]>=20:
       rb=seq[pos-1-rs:pos+n-rs];key=(pos-1,rb,rb[0]);indels[key].add((r.get_tag('RG') if r.has_tag('RG') else '',r.query_name));indread[key]+=1;indstrands[key].add(r.is_reverse)
      pos+=n
     elif op==3:pos+=n
     elif op==4:q+=n
  for (p,ref,alt),frags in indels.items():
   n=indread[(p,ref,alt)];total=depth.get(p,0)
   if len(frags)>=10 and total>=20 and n/total>=.05:
    cands.append({'gene':gene,'chrom':g['chrom'],'pos1':p+1,'ref':ref,'alt':alt,'screen_type':'CIGAR_indel_not_normalized','screen_alt_reads':n,'screen_anchor_read_depth':total,'screen_fraction':n/total,**annotate(gene,p,ref,alt)})
  ds=list(depth.values());summ.append({'gene':gene,'selected_transcript':g['selected_transcript'],'CDS_plus4_bases':len(ds),'minimum_strict_read_depth':min(ds),'bases_ge20':sum(d>=20 for d in ds),'candidates':sum(v['gene']==gene for v in cands)})
  print('screen',gene,summ[-1],round(time.time()-start,1),flush=True);save('screen.checkpoint.json',{'candidates':cands,'gene_summaries':summ})
source=json.loads((OLD/'variants.vcf-records.json').read_text())
for v in cands:
 matches=[r for r in source if (r['CHROM'],r['POS'],r['REF'],r['ALT'])==(v['chrom'],str(v['pos1']),v['ref'],v['alt'])]
 v['Caris_VCF_exact_matches']=[{'line':r['line'],'filter':r['FILTER'],'protein':r['info'].get('PC'),'gene':r['info'].get('GI'),'clinical':r['info'].get('CI'),'VAF':r['sample_format']['VF']} for r in matches]
save('screen-results.json',{'method':'Strict read base-depth screen across selected CDS+/-4bp (22 selected GENCODE37 coding transcripts; tags in target-annotation.json): MAPQ/BQ>=20, duplicate/secondary/supplementary/QCfail/unmapped excluded. SNV initial threshold>=10alt reads and>=5%; indel threshold>=10unique queryname fragments and>=5% of anchor read depth. Indel representations not normalized at screening. Exact source-VCF match performed. Candidate generation only, no somatic or clinical significance claims.','annotation_limit':'Single-variant edit of selected GENCODE37 CDS (21 MANE Select plus NTRK3 APPRIS principal_1), not fully phased patient transcript. Stop-loss outside CDS not assessed; near-exon-boundary effects not comprehensively annotated. In-frame indels are not assigned precise HGVS by this screen.','gene_summaries':summ,'candidates':cands,'elapsed_seconds':time.time()-start})
fields=['gene','chrom','pos1','ref','alt','screen_type','screen_alt_reads','screen_anchor_read_depth','screen_fraction','consequence','protein','coding','transcript','Caris_VCF_exact_matches']
with (B/'screen-candidates.tsv').open('w') as f:
 w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows({k:(json.dumps(v[k]) if k=='Caris_VCF_exact_matches' else v.get(k,'')) for k in fields} for v in cands)

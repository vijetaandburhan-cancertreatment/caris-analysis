"""Dense selected-CDS depth for workbook genes and focal gene neighbors.
No copy-number classification from depth; published targets are unavailable.
"""
import pathlib,json,gzip,re,collections,csv,time,hashlib
import numpy as np,pysam,py2bit
OUT=pathlib.Path(__file__).resolve().parent;ROOT=OUT.parents[2]
SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
ann=json.loads((OUT/'annotation.json').read_text());wb=json.loads((ROOT/'work/oct1-analysis/variants.xlsx-records.json').read_text())
strata={};cnd={}
for r in wb:
 gene=r['Biomarker'].split(' (')[0]
 if r['Test'] in ['Exome Panel - Hybrid - Clinical Genes','Exome Panel - Hybrid - Additional Genes']:strata[gene]='clinical' if 'Clinical' in r['Test'] else 'additional'
 if 'CND' in r['Test']:cnd[gene]={'source_row':r['source_row'],'call':r['Test Result'],'CNA_Value':r['CNA Value']}
selected=set(strata)|set(cnd)
for focal in ['BAP1','RASA1','MTAP','CDKN2A','CDKN2B']:
 g=ann['genes'][focal]
 for gene,a in ann['genes'].items():
  if a['chrom']==g['chrom'] and a['start0']<g['end0']+2000000 and a['end0']>g['start0']-2000000:selected.add(gene)
tx={}
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
 for ln in f:
  if ln.startswith('#'):continue
  a=ln.strip().split('\t')
  if a[2] not in ['transcript','CDS']:continue
  gn=re.search('gene_name "([^"]+)"',a[8]);gn=gn.group(1) if gn else None
  if gn not in selected or gn not in ann['genes'] or a[0]!=ann['genes'][gn]['chrom']:continue
  t=re.search('transcript_id "([^"]+)"',a[8]).group(1);d=tx.setdefault(t,{'gene':gn,'chrom':a[0],'strand':a[6],'transcript':t,'tags':set(),'CDS':[]});d['tags'].update(re.findall('tag "([^"]+)"',a[8]))
  if a[2]=='CDS':d['CDS'].append((int(a[3])-1,int(a[4])))
def rank(t):
 ap=[int(z.rsplit('_',1)[1]) for z in t['tags'] if z.startswith('appris_principal_')]
 return ('MANE_Select'in t['tags'],-min(ap) if ap else -99,'basic' in t['tags'],sum(e-s for s,e in t['CDS']))
bygene=collections.defaultdict(list)
for t in tx.values():
 if t['CDS']:bygene[t['gene']].append(t)
chosen={gene:max(ts,key=rank) for gene,ts in bygene.items()}
for t in chosen.values():t['tags']=sorted(t['tags']);t['CDS']=sorted(t['CDS'])
(OUT/'selected-coding-transcripts.json').write_text(json.dumps(chosen,indent=2)+'\n')
tb=py2bit.open(str(OUT/'public/hg38.2bit'),storeMasked=True);bam=pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai'),threads=1)
def keep(r):return not r.flag&(4|8|256|512|1024|2048) and r.is_proper_pair and r.mapping_quality>=30
def write(p,rows):
 if rows:
  with p.open('w') as f:w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
exons=[];genes=[];started=time.time()
for i,(gene,t) in enumerate(sorted(chosen.items())):
 dd=[];seqall=''
 for j,(s,e) in enumerate(t['CDS']):
  depth=np.sum(np.asarray(bam.count_coverage(t['chrom'],s,e,quality_threshold=25,read_callback=keep)),axis=0);seq=tb.sequence(t['chrom'],s,e);dd.extend(depth.tolist());seqall+=seq
  upper=seq.upper();acgt=sum(upper.count(z) for z in 'ACGT');gc=(upper.count('G')+upper.count('C'))/acgt if acgt else None;rep=sum(z.islower() for z in seq)/len(seq)
  exons.append({'gene':gene,'stratum':strata.get(gene,'not_in_small_variant_workbook'),'chrom':t['chrom'],'start0':s,'end0':e,'CDS_interval_index':j,'transcript':t['transcript'],'length':e-s,'GC':gc,'repeat_fraction':rep,'mean_read_depth':float(depth.mean()),'median_read_depth':float(np.median(depth)),'min_read_depth':int(depth.min()),'max_read_depth':int(depth.max()),'bases_ge20':int((depth>=20).sum()),'bases_ge100':int((depth>=100).sum()),'CND_call':cnd.get(gene,{}).get('call','not_tested')})
 upper=seqall.upper();acgt=sum(upper.count(z) for z in 'ACGT')
 genes.append({'gene':gene,'stratum':strata.get(gene,'not_in_small_variant_workbook'),'chrom':t['chrom'],'gene_start0':ann['genes'][gene]['start0'],'gene_end0':ann['genes'][gene]['end0'],'transcript':t['transcript'],'MANE_Select':int('MANE_Select' in t['tags']),'CDS_intervals':len(t['CDS']),'CDS_bases':len(dd),'GC':(upper.count('G')+upper.count('C'))/acgt if acgt else None,'repeat_fraction':sum(z.islower() for z in seqall)/len(seqall),'mean_read_depth':float(np.mean(dd)),'median_read_depth':float(np.median(dd)),'bases_ge20':sum(x>=20 for x in dd),'CND_call':cnd.get(gene,{}).get('call','not_tested'),'CND_CNA_value':cnd.get(gene,{}).get('CNA_Value','')})
 if i%50==0:print('genes',i+1,'of',len(chosen),'seconds',round(time.time()-started,1),flush=True)
write(OUT/'selected-CDS-interval-depth.tsv',exons);write(OUT/'selected-CDS-gene-depth.tsv',genes)
(OUT/'coding-depth-method.json').write_text(json.dumps({'elapsed_seconds':time.time()-started,'gene_count':len(genes),'CDS_interval_count':len(exons),'selection':'Workbook small-variant and CND genes plus all GENCODE37 protein-coding genes within2Mb ofBAP1/RASA1/MTAP/CDKN2A/B. TranscriptpriorityMANE,APPRISprincipal,basic,longestCDS. PublicCDS is not assaytargetBED.','count':'pysam count_coverage ACGT base depth; MAPQ>=30 BQ>=25 properpaired exclude unmapped/mate-unmapped/duplicate/QCfail/secondary/supplementary. Overlappingmates count separately.','reference':'UCSC hg38.2bit, not exact private customized FASTA','strata':'Workbook test names clinical/additional are provenance descriptors, not proven equal-bait-efficiency capture strata.','CND':'Originalanalyticalworkbooklabelsrecorded for comparison, not independently reclassified. CNA_Value is not assumed calibrated integer tumorcopy number.','limitations':['No normal/captureBED/calibrated reference: adjusted depths must not be called absoluteCN/homozygousloss.','BAP1 small deletion itself removesalignedbases from11bp; that localbasecoverage effect is not exoncopyloss.','Alternate/noncanonical transcripts may be differently captured.']},indent=2)+'\n')
print('complete',len(genes),len(exons),round(time.time()-started,1),flush=True)

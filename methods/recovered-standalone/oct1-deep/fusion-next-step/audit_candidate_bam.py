"""Check selected weak fusion read evidence against original RNA BAM, no calling."""
import collections,gzip,json,pathlib,re
import pysam
r=pathlib.Path(__file__).resolve().parent
p=r/'full/pizzly.insert235.json';d=json.loads(p.read_text());e=json.loads((r/'full/pizzly.insert235.junction-window-validation.json').read_text())
selected={(x['geneA'],x['geneB'])for x in e['summary'] if x['junction40_sequence_pairs']>=2}
gs=[g for g in d['genes'] if (g['geneA']['name'],g['geneB']['name'])in selected]
gids={g[x]['id']for g in gs for x in ['geneA','geneB']};tids={t[x]['id']for g in gs for t in g['transcripts'] for x in ['transcriptA','transcriptB']}
genes={};tx=collections.defaultdict(lambda:{'exon':[],'CDS':[]})
with gzip.open(r.parent/'genomics/gencode.v37.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  v=line.rstrip('\n').split('\t')
  if v[2]not in ['gene','transcript','exon','CDS']:continue
  a=dict(re.findall(r'(\w+) "([^"]*)"',v[8]))
  if v[2]=='gene' and a.get('gene_id')in gids:genes[a['gene_id']]={'chrom':v[0],'start0':int(v[3])-1,'end0':int(v[4]),'strand':v[6],'name':a.get('gene_name')}
  tid=a.get('transcript_id')
  if tid not in tids:continue
  q=tx[tid];q.update(chrom=v[0],strand=v[6],gene_name=a.get('gene_name'),gene_id=a.get('gene_id'))
  if v[2] in ['exon','CDS']:q[v[2]].append((int(v[3]),int(v[4])))
for q in tx.values():
 for k in ['exon','CDS']:q[k].sort(reverse=q['strand']=='-')
def map_nt(tid,pos):
 q=tx[tid];off=0
 for i,(a,b) in enumerate(q['exon']):
  n=b-a+1
  if off<=pos<off+n:
   genomic=a+(pos-off) if q['strand']=='+' else b-(pos-off);cdsoff=0;c=None
   for ca,cb in q['CDS']:
    if ca<=genomic<=cb:c=cdsoff+(genomic-ca if q['strand']=='+' else cb-genomic)+1;break
    cdsoff+=cb-ca+1
   return {'transcript_pos0':pos,'chrom':q['chrom'],'genomic1':genomic,'strand':q['strand'],'exon_number_in_selected_transcript':i+1,'offset_within_exon0':pos-off,'distance_from_exon_end':off+n-1-pos,'coding_position1_if_in_CDS':c}
  off+=n
 return {'error':'position not in annotated exons','transcript_pos0':pos}
needed={x['read1']['name'] for g in gs for x in g['readpairs']};evidence=collections.defaultdict(list);seen=set()
bamfile='/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/RNA_TN26-279853.bam'
index=r.parent.parent/'oct1-analysis/RNA_TN26-279853.bam.bai'
with pysam.AlignmentFile(bamfile,'rb',index_filename=str(index)) as bam:
 for q in genes.values():
  for a in bam.fetch(q['chrom'],q['start0'],q['end0']):
   if a.query_name not in needed:continue
   key=(a.query_name,a.flag,a.reference_id,a.reference_start,a.cigarstring)
   if key in seen:continue
   seen.add(key);fq=a.get_forward_qualities()
   evidence[a.query_name].append({'flag':a.flag,'mate':1 if a.is_read1 else 2,'reference':a.reference_name,'start1':a.reference_start+1,'end1':a.reference_end,'mapq':a.mapping_quality,'cigar':a.cigarstring,'secondary':a.is_secondary,'supplementary':a.is_supplementary,'duplicate_flag':a.is_duplicate,'qc_failed':a.is_qcfail,'query_length':a.query_length,'forward_sequence':a.get_forward_sequence(),'forward_base_qualities':list(fq) if fq is not None else None,'NM':a.get_tag('NM') if a.has_tag('NM') else None,'STAR_nM':a.get_tag('nM') if a.has_tag('nM') else None,'AS':a.get_tag('AS') if a.has_tag('AS') else None,'NH':a.get_tag('NH') if a.has_tag('NH') else None})
out=[]
for g in gs:
 junctions=[]
 for t in g['transcripts']:
  a,b=t['transcriptA'],t['transcriptB']
  junctions.append({'fasta_record':t['fasta_record'],'transcriptA':a,'transcriptB':b,'last_A':map_nt(a['id'],a['endPos']-1),'first_B':map_nt(b['id'],b['startPos'])})
 names={x['read1']['name']for x in g['readpairs']}
 out.append({'geneA':g['geneA'],'geneB':g['geneB'],'pizzly_splitcount':g['splitcount'],'names':sorted(names),'junction_annotation':junctions,'bam_evidence':{n:evidence.get(n,[])for n in sorted(names)}})
result={'source':str(p),'BAM':bamfile,'selected_by':'at least2 distinct sequence-pairs containing the exact central40nt junction window; this is a review priority, not a validated fusion threshold','notes':['Targeted regional BAM fetch can miss a read whose primary alignment lies elsewhere or is unmapped. Missing here does not negate a raw FASTQ candidate.','No RNA duplicate flags/UMIs. Four or more optical-neighbor read names can be one physical template.','Transcript/CDS mapping does not establish biological fusion, translated protein or functional gene inactivation.'],'candidates':out}
(r/'full/candidate-bam-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps([{'pair':g['geneA']['name']+'--'+g['geneB']['name'],'names':len(g['names']),'names_recovered_in_regional_BAM':sum(bool(v) for v in g['bam_evidence'].values()),'first_junction':g['junction_annotation'][0]} for g in out],indent=2))

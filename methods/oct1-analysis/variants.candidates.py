import csv,json,pathlib,collections
p=pathlib.Path('work/oct1-analysis')
v=json.loads((p/'variants.vcf-records.json').read_text())
x=json.loads((p/'variants.xlsx-records.json').read_text())
t=list(csv.DictReader(open('outputs/caris-raw-data/TN26-279853/RNA_TN26-279853.geneTPM_nodup.csv')))
expr={r['Gene']:r for r in t}
lookup={(r['CHROM'],r['POS'],r['info'].get('GI')):r for r in v}
res=[]
for r in x:
    if r['Test Result'] not in ('variantdetected','Pathogenic Variant'):continue
    q=lookup[(r['NGS Chromosome'],r['NGS Position'],r['Biomarker'])]
    ad=q['sample_format']['AD'].split(',')
    res.append({
      'gene':r['Biomarker'],'protein':r['NGS Protein Change'],'coding':r['NGS Coding Change'],
      'transcript':r['NGS Transcript'],'hg38_chrom':q['CHROM'],'position_1based':int(q['POS']),
      'ref':q['REF'],'alt':q['ALT'],'vcf_line':q['line'],'vcf_filter':q['FILTER'],
      'vcf_clinical_impact':q['info'].get('CI','not annotated'),'xlsx_result':r['Test Result'],
      'consequence':q['info'].get('FC'),'tumor_vaf':float(q['sample_format']['VF']),
      'total_depth':int(q['info']['DP']),'reference_reads':int(ad[0]),'variant_reads':int(ad[1]),
      'strand_allele_counts':q['sample_format'].get('SA'),
      'gene_TPM':float(expr.get(r['Biomarker'],{}).get('TPM','nan')),
      'gene_estimated_reads':float(expr.get(r['Biomarker'],{}).get('NumReads','nan')),
      'xlsx_row':r['source_row'],'xlsx_range':'BD%d:CL%d'%(r['source_row'],r['source_row']),
      'clinical_interpretation':r['NGS Interpretation'],
      'somatic_status':'unknown; no matched normal supplied',
      'mutant_RNA_support':'not yet checked; gene TPM does not establish mutant allele expression',
      'research_priority': 'frameshift neoantigen review' if q['info']['FC']=='Frameshift' else ('truncation review; stop alone does not guarantee novel peptide' if q['info']['FC']=='Nonsense' else 'missense review after germline/population filtering')
    })
(p/'variants.candidates.json').write_text(json.dumps(res,indent=2))
with (p/'variants.candidates.tsv').open('w') as f:
    wr=csv.DictWriter(f,res[0].keys(),delimiter='\t');wr.writeheader();wr.writerows(res)
with (p/'variants.candidate-loci.bed').open('w') as f:
    for r in res:
        f.write('%s\t%s\t%s\t%s_%s\n'%(r['hg38_chrom'],max(0,r['position_1based']-101),r['position_1based']+max(100,len(r['ref'])),r['gene'],r['protein']))
selected=['MSLN','CALB2','WT1','PDPN','GATA3','UPK1A','UPK1B','UPK2','UPK3A','UPK3B','EPCAM','CLDN4','NKX2-1','NAPSA','KRT5','KRT6A','KRT7','KRT20','TP63','PAX8','BAP1','NF2','SETD2','MTAP','CDKN2A','CDKN2B','APC','RASA1','PRAME','MAGEA4','MAGEA8','CTAG1B','CTAG2','B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2','TAPBP','PSMB8','PSMB9','CD274','PDCD1','CD3D','CD3E','CD8A','CD8B','IFNG','GZMB','PRF1','CXCL9','CXCL10']
e=[{'gene':g,**expr.get(g,{'TPM':None,'NumReads':None}),'source_csv_row':next((i+2 for i,r in enumerate(t) if r['Gene']==g),None)} for g in selected]
(p/'variants.selected-expression.json').write_text(json.dumps(e,indent=2))
with (p/'variants.selected-expression.tsv').open('w') as f:
    wr=csv.DictWriter(f,['gene','Gene','TPM','NumReads','source_csv_row'],delimiter='\t');wr.writeheader();wr.writerows(e)
print(json.dumps(res,indent=2))

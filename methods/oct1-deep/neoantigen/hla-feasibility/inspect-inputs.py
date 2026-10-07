"""Read-only metadata and supplied expression inspection; no extraction."""
from pathlib import Path
import csv,json,platform,shutil,importlib.metadata,time
import pysam
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
SRC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
OLD=ROOT/'work/oct1-analysis'
bams=[]
for k in ['RNA','DNA']:
    with pysam.AlignmentFile(str(SRC/f'{k}_TN26-279853.bam'),'rb',index_filename=str(OLD/f'{k}_TN26-279853.bam.bai')) as b:
        h=b.header.to_dict();a=[z for z in b.get_index_statistics() if z.contig=='chr6' or 'HLA' in z.contig or z.contig.startswith('chr6_')]
        bams.append({'kind':k,'HD':h['HD'],'PG':[{j:z.get(j) for j in ['PN','VN','ID']} for z in h.get('PG',[])],'index_relevant_counts':[{'reference':z.contig,'mapped':z.mapped,'unmapped':z.unmapped} for z in a],'no_coordinate_count':b.nocoordinate})
expr=[]
with (SRC/'RNA_TN26-279853.geneTPM_nodup.csv').open() as f:
    for z in csv.DictReader(f):
        if z['Gene'].startswith('HLA-') or z['Gene'] in ['B2M','CIITA']:
            expr.append({'gene':z['Gene'],'supplied_gene_TPM':float(z['TPM']),'supplied_estimated_NumReads':float(z['NumReads'])})
hla=[]
for z in json.loads((OLD/'variants.xlsx-records.json').read_text()):
    if z['Biomarker'] in ['HLA-A','HLA-B','HLA-C']:
        hla.append({'gene':z['Biomarker'],'source_workbook_row':z['source_row'],'allele1':z['NGS Protein Change'],'allele2':z['NGS Coding Change'],'method':z['Technology'],'tumor_nuclei_percent':z['HE Percent Tumor Nuclei'],'non_carcinoma_nuclei_percent':z['HE Percent Non Carcinoma Nuclei'],'specimen_type':z['Specimen Type']})
packages={}
for p in ['pysam','biopython','numpy','pandas','scipy','pyarrow']:
    try:packages[p]=importlib.metadata.version(p)
    except importlib.metadata.PackageNotFoundError:packages[p]=None
result={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'architecture':platform.machine(),'free_disk_bytes':shutil.disk_usage(SRC).free,'packages':packages,'PATH_executables':{p:shutil.which(p) for p in ['kallisto','samtools','bedtools','pigz','git-lfs','conda','micromamba']},'bam_index_metadata':bams,'supplied_Caris_classI':hla,'supplied_expression':expr,'limits':['No HLA genotyping, allele-specific quantification, reference installation or read extraction performed.','BAM index counts include aligned records rather than deduplicated molecules and are for sizing only.','Supplied NumReads may be fractional quantification estimates, not unique allele-informative fragments.','Gene expression is bulk tissue; HLA expression can arise from nonmalignant cells.']}
(OUT/'input-feasibility.json').write_text(json.dumps(result,indent=2)+'\n')
with (OUT/'supplied-HLA-gene-expression.tsv').open('w') as f:
    w=csv.DictWriter(f,expr[0].keys(),delimiter='\t');w.writeheader();w.writerows(expr)
print(json.dumps({k:v for k,v in result.items() if k!='supplied_expression'},indent=2))

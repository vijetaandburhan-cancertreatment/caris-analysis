"""Internal DNA/RNA allele consistency, not a validated identity assay."""
from pathlib import Path
from collections import Counter,defaultdict
import json,csv,time
import pysam
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
SRC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
records=json.loads((ROOT/'work/oct1-analysis/variants.vcf-records.json').read_text())
candidates=[r for r in records if len(r['REF'])==len(r['ALT'])==1 and r['ID'].startswith('rs') and set(r['FILTER'].split(';'))<={'rs','Benign','.'} and float(r['sample_format']['VF'])>=.98]
def count(bam,r):
    pos=int(r['POS'])-1;observed=defaultdict(set);excluded=Counter()
    for read in bam.fetch(r['CHROM'],pos,pos+1):
        if read.flag&(4|256|512|1024|2048) or read.mapping_quality<30 or read.query_qualities is None:continue
        q=next((q for q,p in read.get_aligned_pairs(matches_only=True) if p==pos),None)
        if q is None or read.query_qualities[q]<30:continue
        if min(q-read.query_alignment_start,read.query_alignment_end-1-q)<5:continue
        if read.has_tag('NH') and read.get_tag('NH')!=1:continue
        key=(read.get_tag('RG') if read.has_tag('RG') else '',read.query_name)
        observed[key].add(read.query_sequence[q])
    calls=Counter(next(iter(a)) if len(a)==1 else 'discordant' for a in observed.values())
    denominator=sum(calls[a] for a in 'ACGT')
    return {'counts':dict(calls),'assessable_names':denominator,'specified_alt_fraction':calls[r['ALT']]/denominator if denominator else None}
rows=[];started=time.time()
handles={kind:pysam.AlignmentFile(str(SRC/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai')) for kind in ['DNA','RNA']}
for r in candidates:
    row={'chrom':r['CHROM'],'pos1':int(r['POS']),'public_id':r['ID'],'ref':r['REF'],'alt':r['ALT'],'gene':r['info'].get('GI'),'Caris_filter':r['FILTER'],'Caris_VAF':r['sample_format']['VF']}
    row.update({kind:count(bam,r) for kind,bam in handles.items()});rows.append(row)
for b in handles.values():b.close()
eligible=[r for r in rows if r['DNA']['assessable_names']>=20 and r['RNA']['assessable_names']>=20 and r['DNA']['specified_alt_fraction']>=.98]
discordant=[r for r in eligible if r['RNA']['specified_alt_fraction']<.90]
summary={'tested_export_selected_markers':len(rows),'DNA_high_alt_and_both_depth20_markers':len(eligible),'RNA_same_allele_ge90pct':sum(r['RNA']['specified_alt_fraction']>=.90 for r in eligible),'RNA_same_allele_ge98pct':sum(r['RNA']['specified_alt_fraction']>=.98 for r in eligible),'discordant_below90pct':discordant,'genes_represented':len(set(r['gene'] for r in eligible)),'chromosomes_represented':sorted(set(r['chrom'] for r in eligible)),'elapsed_seconds':time.time()-started}
result={'scope':'Internal consistency of export-selected dbSNP alleles with reported DNA VAF >=98%, excluding export low-quality filters. Recount uses primary nonduplicate MAPQ30/BQ30, bases >=5 from aligned ends; RNA unique mappings when NH present; mates collapsed and discordants excluded. Not an ancestry/identity test, not independently selected fingerprint loci, not matched-normal analysis, not contamination estimate. Selection from this same DNA VCF creates ascertainment bias; no chance-match probability calculated.','summary':summary,'rows':rows}
(OUT/'findings.json').write_text(json.dumps(result,indent=2)+'\n')
with (OUT/'allele-consistency.tsv').open('w') as f:
    fields=['chrom','pos1','public_id','ref','alt','gene','DNA_names','DNA_alt_fraction','RNA_names','RNA_alt_fraction']
    w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader()
    for r in rows:w.writerow({**{k:r[k] for k in fields[:6]},'DNA_names':r['DNA']['assessable_names'],'DNA_alt_fraction':r['DNA']['specified_alt_fraction'],'RNA_names':r['RNA']['assessable_names'],'RNA_alt_fraction':r['RNA']['specified_alt_fraction']})
print(json.dumps(summary,indent=2),flush=True)

"""Conservative bulk TCR summary and independent exact CDR3 read recount.

No cell barcodes, chain pairing, antigen specificity, clonality score or tumor
reactivity is inferred. Local patient sequences are written only to local files.
"""
from pathlib import Path
from collections import Counter,defaultdict
from Bio.Seq import Seq
import csv,hashlib,itertools,json,re
import pysam
OUT=Path(__file__).resolve().parent
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/trust4-tools')
PAT=BASE/'patient';prefix='RNA_TN26-279853'
run=json.loads((OUT/'patient-rna-result.json').read_text());assert run['status']=='completed'
report=PAT/(prefix+'_report.tsv');cdr3=PAT/(prefix+'_cdr3.out');airr=PAT/(prefix+'_airr.tsv')
rows=list(csv.DictReader(report.open(),delimiter='\t'))
scores=defaultdict(list)
for line in cdr3.read_text().splitlines():
    a=line.split()
    if len(a)>=13:scores[(a[0],a[8])].append({'CDR3_score':float(a[9]),'assigned_fragment_count':float(a[10]),'CDR3_germline_similarity':float(a[11]),'complete_vdj_assembly':a[12]})
def chain(row):
    for key in ('C','J','V'):
        s=row[key]
        if s.startswith(('TRA','TRB','TRD','TRG','IGH','IGK','IGL')):return s[:3]
    return 'unknown'
selected=[];allchains=Counter()
for r in rows:
    ch=chain(r);allchains[ch]+=1
    if ch not in ('TRA','TRB','TRD','TRG'):continue
    nt=r['CDR3nt'];valid=bool(re.fullmatch('[ACGT]+',nt)) and len(nt)%3==0
    aa=str(Seq(nt).translate()) if valid else None
    sc=scores.get((r['cid'],nt),[])
    selected.append({'research_candidate_id':f'TCR{len(selected)+1:04d}','chain':ch,'CDR3_nt':nt,'CDR3_aa_reported':r['CDR3aa'],'CDR3_aa_independently_translated':aa,'V':r['V'],'D':r['D'],'J':r['J'],'C':r['C'],'consensus_id':r['cid'],'tool_read_count':float(r['#count']),'tool_report_frequency':float(r['frequency']),'tool_complete_vdj':r['cid_full_length'],'cdr3_details':sc,'max_CDR3_score':max([x['CDR3_score'] for x in sc],default=None),'unambiguous_in_frame_stop_free_CDR3':bool(valid and '*' not in aa and aa==r['CDR3aa'])})
targets=defaultdict(list)
for i,r in enumerate(selected):
    if not re.fullmatch('[ACGT]+',r['CDR3_nt']):continue
    nt=r['CDR3_nt']
    for orientation,pat in [('forward',nt),('reverse',str(Seq(nt).reverse_complement()))]:
        targets[pat[:15]].append((i,orientation,pat))
support=[{'all_exact_fragment_names':set(),'Q20_fragment_names':set(),'Q30_fragment_names':set(),'Q20_5base_margin_fragment_names':set(),'read_sequence_sha256':set(),'Q20_read_sequence_sha256':set(),'exact_reads':0,'Q20_reads':0,'forward_match_reads':0,'reverse_match_reads':0} for _ in selected]
fastqs=[PAT/(prefix+'_toassemble_1.fq'),PAT/(prefix+'_toassemble_2.fq')]
fqcounts=Counter();names=set();name_digest=hashlib.sha256();paired_hash=[hashlib.sha256(),hashlib.sha256()]
expected_source={};expected_candidate_indices=defaultdict(set)
def norm(name):return re.sub('/[12]$','',name)
with pysam.FastxFile(str(fastqs[0])) as one,pysam.FastxFile(str(fastqs[1])) as two:
    for a,b in itertools.zip_longest(one,two):
        assert a is not None and b is not None and norm(a.name)==norm(b.name)
        qname=norm(a.name);fqcounts['complete_paired_records']+=1
        if qname in names:fqcounts['duplicate_pair_names']+=1
        names.add(qname);name_digest.update((qname+'\n').encode())
        for mate,r in enumerate((a,b)):
            assert r.quality is not None and len(r.sequence)==len(r.quality)
            seq=r.sequence.upper();qual=r.quality
            paired_hash[mate].update(('@'+r.name+'\n'+seq+'\n+\n'+qual+'\n').encode())
            fqcounts['reads']+=1;fqcounts['bases']+=len(seq)
            found={}
            for pos in range(max(0,len(seq)-14)):
                for i,orientation,pat in targets.get(seq[pos:pos+15],[]):
                    if seq.startswith(pat,pos):
                        end=pos+len(pat);q=min(ord(x)-33 for x in qual[pos:end])
                        margin=min(pos,len(seq)-end)
                        if i not in found or q>found[i][0]:found[i]=(q,margin,orientation)
            fingerprint=hashlib.sha256(seq.encode()).hexdigest()
            for i,(q,margin,orientation) in found.items():
                s=support[i];s['all_exact_fragment_names'].add(qname);s['exact_reads']+=1;s['read_sequence_sha256'].add(fingerprint)
                s[orientation+'_match_reads']+=1
                if q>=20:
                    s['Q20_fragment_names'].add(qname);s['Q20_reads']+=1;s['Q20_read_sequence_sha256'].add(fingerprint)
                    if margin>=5:s['Q20_5base_margin_fragment_names'].add(qname)
                    expected_source[(qname,mate+1)]=(seq,qual)
                    expected_candidate_indices[(qname,mate+1)].add(i)
                if q>=30:s['Q30_fragment_names'].add(qname)
assert fqcounts['duplicate_pair_names']==0
source_counts=Counter();verified_mates=set();mismatched=[]
src=BASE.parent/'TN26-279853/RNA_TN26-279853.bam'
idx=OUT.parents[3]/'work/oct1-analysis/RNA_TN26-279853.bam.bai'
with pysam.AlignmentFile(str(src),'rb',index_filename=str(idx)) as bam:
    source_counts['header_read_groups']=len(bam.header.to_dict().get('RG',[]))
    assert source_counts['header_read_groups']==1, 'Need RG-aware source verification for multiple groups'
    for region in ['chr7','chr14','*']:
        for r in bam.fetch(region):
            if r.flag&(256|2048):continue
            k=(r.query_name,1 if r.is_read1 else 2 if r.is_read2 else 0)
            if k not in expected_source:continue
            source_counts['matching_primary_records_found']+=1
            if r.query_qualities is None:source_counts['missing_source_quality']+=1;continue
            seq=r.get_forward_sequence();qual=''.join(chr(q+33) for q in r.get_forward_qualities())
            if (seq,qual)==expected_source[k]:verified_mates.add(k)
            else:
                source_counts['sequence_mismatches']+=int(seq!=expected_source[k][0]);source_counts['quality_mismatches']+=int(qual!=expected_source[k][1]);mismatched.append(k)
source_counts['Q20_supporting_mates_expected']=len(expected_source)
source_counts['Q20_supporting_mates_exactly_verified_in_original_BAM']=len(verified_mates)
source_counts['Q20_supporting_mates_not_found_in_chr7_chr14_unmapped_scan']=len(set(expected_source)-verified_mates-set(mismatched))
assert not mismatched, 'Unexpected source sequence/quality mismatch; stop rather than claim Q20 evidence'
verified_fragments=[set() for _ in selected]
for k in verified_mates:
    for i in expected_candidate_indices[k]:verified_fragments[i].add(k[0])
for r,s in zip(selected,support):
    for k,v in s.items():r['independent_'+k.replace('_names','_count').replace('_sha256','_unique_count')]=len(v) if isinstance(v,set) else v
    r['conservative_research_shortlist']=bool(r['chain'] in ('TRA','TRB') and r['unambiguous_in_frame_stop_free_CDR3'] and r['V']!='.' and r['J']!='.' and r['max_CDR3_score'] is not None and r['max_CDR3_score']>=0.5 and len(s['Q20_fragment_names'])>=2 and len(s['Q20_read_sequence_sha256'])>=2)
    # Hash names, never publish read names in findings.
    r['exact_fragment_name_set_sha256']=hashlib.sha256('\n'.join(sorted(s['all_exact_fragment_names'])).encode()).hexdigest()
for i,r in enumerate(selected):
    r['original_BAM_verified_Q20_fragment_count']=len(verified_fragments[i])
    r['conservative_research_shortlist']=r['conservative_research_shortlist'] and len(verified_fragments[i])>=2
selected.sort(key=lambda r:(not r['conservative_research_shortlist'],-r['independent_Q20_fragment_count'],-r['tool_read_count']))
fields=[k for k in selected[0] if k!='cdr3_details'] if selected else ['research_candidate_id','chain']
for filename,data in [('all-reconstructed-TCR.tsv',selected),('read-supported-TRA-TRB.tsv',[r for r in selected if r['conservative_research_shortlist']])]:
    with (OUT/filename).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(data)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
summary={'status':'exploratory_complete','tool':'TRUST4 v1.1.11-r641','source_commit':'a3fedd4aa0c1ad4da815427d82e0ceb69ce9c3a0','reference_files':'Unchanged bundled hg38_bcrtcr.fa and human_IMGT+C.fa, exact hashes in public-source-manifest.json','run':run,'official_control':json.loads((OUT/'public-control-comparison.json').read_text()),'all_reported_chain_rows':dict(allchains),'total_reported_rows':len(rows),'TCR_rows':len(selected),'conservative_TRA_TRB_rows':sum(r['conservative_research_shortlist'] for r in selected),'conservative_chain_counts':dict(Counter(r['chain'] for r in selected if r['conservative_research_shortlist'])),'unique_conservative_chain_CDR3_nt':len({(r['chain'],r['CDR3_nt']) for r in selected if r['conservative_research_shortlist']}),'extracted_FASTQ_check':{'counts':dict(fqcounts),'unique_pair_names':len(names),'ordered_pair_name_sha256':name_digest.hexdigest(),'FASTQ_content_sha256':[h.hexdigest() for h in paired_hash]},'candidates':selected,'file_hashes':[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in [report,cdr3,airr,*fastqs] if p.exists()], 'definitions':{'CDR3_productivity_proxy':'CDR3 nucleotide sequence contains only A/C/G/T, length divisible by three, translates exactly to reported amino acids without stop codon. This is a CDR3-local in-frame/stop-free proxy, not demonstration of a functional full receptor.','direct_read_support':'Full exact CDR3 nucleotide substring (either orientation) in extracted patient RNA reads; minimum base Q20 across the whole CDR3 for stringent count; collapse mates by original query name. No MAPQ filter because rearranged reads may be unmapped. Distinct full read-sequence hashes add position/sequence diversity context, not independent molecules.','conservative_shortlist':'TRA or TRB, V and J annotated, local in-frame/stop-free CDR3, CDR3_score>=0.5 (excludes imputed0.01 and partial0), >=2 exact Q20 paired-name fragments and >=2 different Q20-supporting full-read sequences. Complete VDJ assembly is reported separately, not required for a CDR3 comparison resource.'},'limitations':['No cell barcodes or paired-chain inference: TRA and TRB rows must not be combined into receptors.','No antigen specificity, tumor reactivity, patient-cell identity, validated repertoire clonality/diversity, or treatment-response inference.','Bulk FFPE brain tumor includes nonmalignant cells; recovered RNA may be from infiltrating or circulating lymphocytes.','Counts and frequencies are sequencing/assembly measures, not cell proportions; no UMI/PCR deduplication and incomplete sampling/FFPE damage can bias them.','Gene/allele assignments may be ambiguous. Exact public-control CDR3 identity/count agreement is execution validation; the BCR-only control is not a TCR sensitivity benchmark.','Some chains may have only complete CDR3 rather than complete VDJ; do not use partial reconstructed chains directly as engineered TCR constructs.','A negative result would not establish absence of T cells or tumor-reactive TILs.','Optional value is comparing exact CDR3 nucleotide/amino-acid and V/J assignments with future validated paired single-cell TCR data; a match alone would not establish specificity.'],'sources':[{'title':'Official TRUST4 source and workflow','url':'https://github.com/liulab-dfci/TRUST4/tree/a3fedd4aa0c1ad4da815427d82e0ceb69ce9c3a0'},{'title':'TRUST4 primary methods study','url':'https://www.nature.com/articles/s41592-021-01142-2'}]}
summary['source_BAM_verification']={'regions':['chr7','chr14','no-coordinate/unmapped'],'counts':dict(source_counts),'method':'Every Q20-supporting read found in these indexed regions was compared byte-for-byte for original-orientation sequence and quality against the source BAM primary record. Additional conservative rule requires >=2 source-verified supporting fragment names. Unfound mates are explicit and cannot count toward that rule.'}
summary['definitions']['conservative_shortlist']+=' Also requires >=2 direct Q20 supporting fragment names with exact sequence/quality verification back to the original BAM.'
(OUT/'findings.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:summary[k] for k in ['status','all_reported_chain_rows','total_reported_rows','TCR_rows','conservative_TRA_TRB_rows','conservative_chain_counts','unique_conservative_chain_CDR3_nt','extracted_FASTQ_check']},indent=2))

"""Bounded BAP1 peptide nucleotide/local-splicing validation, not clinical calling.

No sequence leaves this machine. Original files opened read-only. Query-name
fragment collapse joins mates but is NOT UMI/PCR molecule deduplication.
"""
from pathlib import Path
import csv, collections, gzip, hashlib, json, re, statistics, time
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[3]
SRC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
OLD=ROOT/'work/oct1-analysis'
GEN=ROOT/'work/oct1-deep/genomics'
CH='chr3'; MASK=4|256|512|1024|2048
def write(name,data): (OUT/name).write_text(json.dumps(data,indent=2)+'\n')
def tsv(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,keys,delimiter='\t');w.writeheader();w.writerows(rows)
def ok(r):return not r.flag&MASK and r.mapping_quality>=20 and r.query_qualities is not None
def key(r):return (r.get_tag('RG') if r.has_tag('RG') else '',r.query_name)
def collapse(d):return dict(collections.Counter(next(iter(v)) if len(v)==1 else 'discordant' for v in d.values()))
def layout(r):
    p=r.reference_start;q=0;b={};d=[];ins=[];j=[]
    for op,n in r.cigartuples or []:
        if op in (0,7,8):b.update((p+i,q+i) for i in range(n));p+=n;q+=n
        elif op==1:ins.append((p,q,n));q+=n
        elif op==2:d.append((p,n));p+=n
        elif op==3:j.append((p,p+n));p+=n
        elif op==4:q+=n
        elif op in (5,6):pass
        else:raise ValueError(op)
    return b,d,ins,j

annot=json.loads((GEN/'target-annotation.json').read_text())['genes']['BAP1']
record_path=OLD/'transcripts/NM_004656.3.gb'
record=SeqIO.read(record_path,'genbank')
f=next(f for f in record.features if f.type=='CDS')
cds=str(f.extract(record.seq));protein=str(Seq(cds).translate(cds=True))
assert protein==f.qualifiers['translation'][0]
cpos=[p for st,en in sorted(annot['selected_CDS'],reverse=True) for p in range(en-1,st-1,-1)]
g2c={p:i+1 for i,p in enumerate(cpos)}
assert len(cpos)==len(cds)-3
# Independent exact public hg38 target-exon sequence versus RefSeq projection.
rw=next(r for r in json.loads((OLD/'variant-read-support.hg38-reference-windows.json').read_text()) if r['gene']=='BAP1')
ref=rw['dna'].upper(); rs=rw['start'];exst=52408473;exen=52408606
assert str(Seq(ref[exst-rs:exen-rs]).reverse_complement())==cds[122:255]
delst=52408550;delen=11
edited=ref[:delst-rs]+ref[delst-rs+delen:]
eq=[s for s in range(rs+3,rw['end']-delen-3) if ref[:s-rs]+ref[s-rs+delen:]==edited]
assert delst in eq
mutcds=cds[:167]+cds[178:]
peps=[]
for label,a,b,ma,mb,expected in [
    ('IEERKGLYL',157,194,157,183,'IEERKGLYL'),
    ('EERKGLYL',160,194,160,183,'EERKGLYL'),
    ('extended_context_through_predicted_stop',148,203,148,192,'FKWIEERKGLYLGG*')]:
    st=cpos[b-1];en=cpos[a-1]+1
    wt=ref[st-rs:en-rs];mt=wt[:delst-st]+wt[delst-st+delen:]
    assert str(Seq(mt).reverse_complement())==mutcds[ma-1:mb]
    assert str(Seq(mt).reverse_complement().translate())==expected
    peps.append({'label':label,'ref_c_start':a,'ref_c_end':b,'mutant_c_start':ma,'mutant_c_end':mb,'start0':st,'end0':en,'expected_peptide':expected,'ref_nt':wt,'mutant_nt':mt})
exons=[]
for line in gzip.open(GEN/'gencode.v37.annotation.gtf.gz','rt'):
    if ('transcript_id "'+annot['selected_transcript']+'"') in line and '\texon\t' in line:
        a=line.rstrip().split('\t');exons.append({'number':int(re.search(r'exon_number (\d+)',a[8]).group(1)),'start0':int(a[3])-1,'end0':int(a[4])})
exons.sort(key=lambda x:x['number'])
known_j={ (exons[i+1]['end0'],exons[i]['start0']):f'exon_{i+1}_to_{i+2}' for i in range(len(exons)-1)}
skipj=(exons[4]['end0'],exons[2]['start0'])
assert known_j[(52408606,52409553)]=='exon_3_to_4'
assert known_j[(52408077,52408473)]=='exon_4_to_5'
region=(52407957,52409878) # selected coding exons 1--5, c.1--375
coding=[p for p in cpos[:375]]
expectedbase={p:str(Seq(cds[c-1]).complement()) for p,c in g2c.items() if p in coding}
vcf=[]
for line in (SRC/'DNA_TN26-279853.vcf').open():
    if line.startswith('#'):continue
    a=line.rstrip().split('\t')
    if a[0]==CH and region[0]<=int(a[1])-1<region[1]:
        info=dict(z.split('=',1) if '=' in z else (z,True) for z in a[7].split(';'))
        vcf.append({'chrom':a[0],'pos1':int(a[1]),'ref':a[3],'alt':a[4],'filter':a[6],'coding':info.get('DC'),'protein':info.get('PC')})
write('reference-provenance.json',{'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'RefSeq':record.id,'RefSeq_file':str(record_path),'RefSeq_sha256':hashlib.sha256(record_path.read_bytes()).hexdigest(),'selected_MANE_GENCODE37':annot['selected_transcript'],'GTF_sha256':json.loads((GEN/'target-annotation.json').read_text())['sha256'],'reference_build':'hg38','reference_window':rw,'target_exon':{'number':4,'start0':exst,'end0':exen,'coding_start1':123,'coding_end1':255},'target_exon_RefSeq_public_hg38_exact_match':True,'genomic_equivalent_deletion_starts0':eq,'peptide_windows':peps,'all_selected_exons':exons,'source_VCF_nearby_calls':vcf})

allresults=[];mismatchrows=[];junctionrows=[];perbaserows=[];indelrows=[]
for kind in ['DNA','RNA']:
    with pysam.AlignmentFile(str(SRC/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(OLD/f'{kind}_TN26-279853.bam.bai')) as bam:
        # Preserve filter-denominator detail; reads fetched once across five exons.
        reads=[];bases_bypos=collections.defaultdict(collections.Counter); basefrag=collections.defaultdict(lambda:collections.defaultdict(set))
        df=collections.defaultdict(set); dreads=collections.Counter(); dsame={};ind=collections.defaultdict(set);incov=collections.defaultdict(set)
        for r in bam.fetch(CH,*region):
            if not ok(r):continue
            b,d,i,j=layout(r);k=key(r);call='unassessable'
            # Three matching, Q20 aligned bases either side of an equivalent del.
            for s,n in d:
                if n==delen and s in eq:
                    pp=list(range(s-3,s))+list(range(s+n,s+n+3))
                    if all(p in b and r.query_qualities[b[p]]>=20 and r.query_sequence[b[p]]==ref[p-rs] for p in pp):call='alt'
            dst=min(eq)-3;den=max(eq)+delen+3
            if call!='alt' and all(p in b and r.query_qualities[b[p]]>=20 and r.query_sequence[b[p]]==ref[p-rs] for p in range(dst,den)) and not any(dst<s<den for s,q,n in i):call='ref'
            if call in ('alt','ref'):df[k].add(call);dreads[call]+=1
            reads.append((r,b,d,i,j,call));dsame[id(r)]=call
            for p,q in b.items():
                if p in expectedbase and r.query_qualities[q]>=20:
                    base=r.query_sequence[q];bases_bypos[p][base]+=1;basefrag[p][k].add(base)
            for s,n in d:
                if s in expectedbase or s+n-1 in expectedbase:
                    pp=list(range(s-3,s))+list(range(s+n,s+n+3))
                    if all(p in b and r.query_qualities[b[p]]>=20 for p in pp):ind[(s,'D',n)].add(k)
            for s,q,n in i:
                if s in expectedbase or s-1 in expectedbase:
                    pp=list(range(s-3,s))+list(range(s,s+3))
                    if all(p in b and r.query_qualities[b[p]]>=20 for p in pp) and min(r.query_qualities[q:q+n])>=20:ind[(s,'I',n)].add(k)
        dfrag=collapse(df)
        # High-quality per-position observations and variants above predeclared
        # modest research thresholds. Low fraction counts retained separately.
        candidates=[]
        for p in sorted(expectedbase):
            rc=bases_bypos[p];fc=collapse(basefrag[p]);refb=expectedbase[p]
            den=sum(n for a,n in fc.items() if a in 'ACGT');alts={a:n for a,n in fc.items() if a in 'ACGT' and a!=refb}
            row={'kind':kind,'pos1':p+1,'coding_pos1':g2c[p],'ref_forward':refb,'depth_reads':sum(rc.values()),'depth_fragments':den,'ref_fragments':fc.get(refb,0),'alt_fragments_total':sum(alts.values()),'discordant_mate_fragments':fc.get('discordant',0),'base_read_counts':json.dumps(dict(rc),sort_keys=True),'base_fragment_counts':json.dumps(fc,sort_keys=True)}
            perbaserows.append(row)
            for a,n in alts.items():
                if n>=(5 if kind=='DNA' else 3) and den and n/den>=.01:
                    ondel=collections.Counter();single=collections.Counter();strand=collections.Counter();coords=set()
                    for r,b,d,i,j,dc in reads:
                        if p in b and r.query_qualities[b[p]]>=20 and r.query_sequence[b[p]]==a:
                            strand['reverse' if r.is_reverse else 'forward']+=1;coords.add((r.reference_start,r.reference_end,r.cigarstring,r.is_reverse,r.is_read1))
                            if dc in ('alt','ref'):single[dc]+=1
                    for k,z in basefrag[p].items():
                        if z=={a}:
                            z2=df.get(k,set());ondel[next(iter(z2)) if len(z2)==1 else 'discordant' if z2 else 'unphased']+=1
                    cand={**row,'alt_forward':a,'alt_fragments':n,'alt_fraction':n/den,'linked_deletion_fragment_counts':dict(ondel),'same_read_deletion_counts':dict(single),'alt_read_strands':dict(strand),'distinct_alt_alignment_patterns':len(coords)}
                    candidates.append(cand);mismatchrows.append(cand)
        # Indel observations with low-level support retained; ascertainment
        # is CIGAR-based, so these are not normalized/validated variant calls.
        for (s,op,n),ks in sorted(ind.items()):
            if len(ks)>=3:
                indelrows.append({'kind':kind,'position0':s,'operation':op,'length':n,'query_name_fragments':len(ks),'is_equivalent_known_BAP1_deletion':op=='D' and n==delen and s in eq,'coding_pos1_at_start':g2c.get(s),'phase_to_known_deletion':dict(collections.Counter(next(iter(df[k])) if len(df.get(k,set()))==1 else 'unphased_or_discordant' for k in ks))})
        local=[]
        for win in peps:
            st=win['start0'];en=win['end0'];cnt=collections.Counter();fr=collections.defaultdict(set);strands=collections.Counter();coords=collections.defaultdict(set);clean=collections.defaultdict(set);translated=collections.defaultdict(set);qc=collections.Counter();del_observed=collections.Counter();maps=collections.defaultdict(list)
            for r,b,d,i,j,dc in reads:
                if st not in b or en-1 not in b:continue
                if any(s<en and e>st for s,e in j):qc['skipped_window']+=1;continue
                qa=b[st];qb=b[en-1]+1
                if min(r.query_qualities[qa:qb])<20:qc['low_full_window_BQ']+=1;continue
                obs=r.query_sequence[qa:qb]
                call='exact_deletion_haplotype' if obs==win['mutant_nt'] else 'exact_reference_haplotype' if obs==win['ref_nt'] else 'other_haplotype'
                cnt[call]+=1;k=key(r);fr[k].add(call);strands[call+('_reverse' if r.is_reverse else '_forward')]+=1
                coords[call].add((r.reference_start,r.reference_end,r.cigarstring,r.is_reverse,r.is_read1));maps[call].append(r.mapping_quality)
                if not any(op==4 for op,n in r.cigartuples) and min(qa-r.query_alignment_start,r.query_alignment_end-qb)>=5:clean[call].add(k)
                if any(s in eq and n==delen for s,n in d):
                    del_observed[call]+=1
                    # Only a complete in-frame window with expected deletion
                    # length is translated. No inference from partial overlaps.
                    if len(obs)==len(win['mutant_nt']):translated[k].add(str(Seq(obs).reverse_complement().translate()))
            # Exclude mates disagreeing on nucleotide-haplotype class from the
            # clean subset, and exclude translation-discordant mates explicitly.
            lc={'label':win['label'],'reads':dict(cnt),'fragments':collapse(fr),'read_strands':dict(strands),'distinct_alignment_patterns':{k:len(v) for k,v in coords.items()},'clean_no_softclip_at_least_5bp_from_ends_fragments':{a:sum(fr[k]=={a} for k in v) for a,v in clean.items()},'median_MAPQ':{k:statistics.median(v) for k,v in maps.items()},'window_quality_exclusions':dict(qc),'reads_containing_equivalent_deletion_by_haplotype':dict(del_observed),'deletion_reads_predicted_translation_fragment_counts':collapse(translated)}
            local.append(lc)
        junction_summary=[]
        if kind=='RNA':
            js=collections.defaultdict(lambda:{'reads':0,'f':set(),'coords':set(),'same_read_del':collections.Counter(),'freads':0,'rreads':0,'proper_pair_f':set()})
            for r,b,d,i,j,dc in reads:
                for s,e in j:
                    # Include boundaries adjoining exon 4, and complete skips.
                    if not (s<=exen and e>=exst):continue
                    pp=list(range(s-8,s))+list(range(e,e+8))
                    if not all(p in b and r.query_qualities[b[p]]>=20 for p in pp):continue
                    if [b[p] for p in range(s-8,s)]!=list(range(b[s-8],b[s-8]+8)) or [b[p] for p in range(e,e+8)]!=list(range(b[e],b[e]+8)):continue
                    z=js[(s,e)];z['reads']+=1;z['f'].add(key(r));z['coords'].add((r.reference_start,r.reference_end,r.cigarstring,r.is_reverse,r.is_read1));z['same_read_del'][dc]+=1;z['rreads' if r.is_reverse else 'freads']+=1
                    if r.is_proper_pair:z['proper_pair_f'].add(key(r))
            for (s,e),z in sorted(js.items()):
                phase=collections.Counter();phaseproper=collections.Counter()
                for k in z['f']:
                    calls=df.get(k,set());label=next(iter(calls)) if len(calls)==1 else 'discordant' if calls else 'unphased';phase[label]+=1
                    if k in z['proper_pair_f']:phaseproper[label]+=1
                jr={'kind':'RNA','intron_start0':s,'intron_end0':e,'annotation':known_j.get((s,e),'exon_3_to_5_skip_of_exon_4' if (s,e)==skipj else 'other_junction'),'reads':z['reads'],'fragments':len(z['f']),'distinct_alignment_patterns':len(z['coords']),'forward_reads':z['freads'],'reverse_reads':z['rreads'],'same_read_deletion_state':dict(z['same_read_del']),'same_query_fragment_deletion_state':dict(phase),'proper_pair_query_fragment_deletion_state':dict(phaseproper)}
                junction_summary.append(jr);junctionrows.append(jr)
            upper=js.get((52408606,52409553),{'f':set()})['f'];lower=js.get((52408077,52408473),{'f':set()})['f'];both=upper&lower
            both_summary={'fragments_with_both_exon4_junctions':len(both),'deletion_state':dict(collections.Counter(next(iter(df[k])) if len(df.get(k,set()))==1 else 'unphased_or_discordant' for k in both)),'exon4_skip_junction_fragments':len(js.get(skipj,{'f':set()})['f'])}
        else:both_summary=None
        exon_depth=[]
        for num in range(1,6):
            ex=exons[num-1];ps=[p for p in expectedbase if ex['start0']<=p<ex['end0']];dep=[sum(bases_bypos[p].values()) for p in ps]
            exon_depth.append({'exon':num,'coding_bases':len(ps),'median_read_depth':statistics.median(dep),'min_read_depth':min(dep),'max_read_depth':max(dep)})
        result={'kind':kind,'strict_deletion_reads':dict(dreads),'strict_deletion_fragments':dfrag,'screened_coding_positions':len(expectedbase),'exon_depth':exon_depth,'nearby_SNV_candidates_above_threshold':candidates,'peptide_haplotype_support':local,'junctions':junction_summary,'dual_exon4_junction_support':both_summary}
        allresults.append(result)
        print(kind,'deletion',dfrag,'nearby_SNV_candidates',len(candidates),'local peptide',[(x['label'],x['fragments']) for x in local],flush=True)

limits=['RNA/DNA nucleotide and splice-read evidence only: no translated protein, processing, HLA presentation, T-cell recognition or therapeutic activity measured.','Mutant peptide translation assumes the named RefSeq CDS initiation and reading frame; local read phasing is not full-length transcript phasing.','Query-name+read-group fragments merge mates, not PCR duplicates or independent UMI-tagged molecules. RNA BAM absence of duplicate flags is not proof of deduplication.','No matched normal: local allele support does not establish tumor specificity, clonality, germline status or BAP1 second-hit mechanism.','Nearby screen is bounded to first five MANE coding exons with modest exploratory thresholds; it is not a clinical negative call or de novo complete variant callset.','RNA MAPQ255 is the aligner unique-mapping sentinel, not measured phred255 confidence.','Splice junction counts require 8 aligned contiguous Q20 bases on each side; shorter anchors and low-quality/supplementary evidence are excluded.','RNA abundance/alternate reads do not quantify nonsense-mediated decay; a premature stop prediction is not proof of stability or protein production.']
write('results.json',{'method':{'flags_excluded':MASK,'minimum_MAPQ':20,'minimum_BQ':20,'nearby_SNV_threshold':'At least 5 DNA or 3 RNA non-discordant query-name fragments and >=1% at a position; threshold for review, not validated variant sensitivity.','peptide_haplotype':'Each qualifying individual read must span every peptide nucleotide and the full deletion, Q20 throughout; exact genomic haplotype compared before reverse-complement translation.','splicing':'CIGAR N junctions, contiguous 8bp Q20 anchors on both sides, same strict flags/MAPQ. Phase on same read or RG+query-name separately.','screen_coding_region':'MANE coding exons 1--5, c.1--375; selected RefSeq sequence projected onto GENCODE v37 MANE CDS coordinates, target exon independently matches existing public hg38 window.'},'results':allresults,'limitations':limits})
tsv('nearby-coding-base-evidence.tsv',perbaserows);tsv('nearby-SNV-review-candidates.tsv',mismatchrows);tsv('local-indel-observations.tsv',indelrows);tsv('exon4-splice-junction-evidence.tsv',junctionrows)

#!/usr/bin/env python3
import pathlib,importlib.util,json,array,pysam
W=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('production',W/'audit_target_splicing.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
results=[]
def read(seq,cigar,start=100,flag=0,quality=30):
 r=pysam.AlignedSegment();r.query_name='synthetic';r.query_sequence=seq;r.flag=flag;r.reference_id=0;r.reference_start=start;r.mapping_quality=60;r.cigartuples=cigar;r.query_qualities=array.array('B',[quality]*len(seq));return r
def check(name,observed,expected):
 results.append({'name':name,'observed':observed,'expected':expected,'pass':observed==expected})
genome='GTCGATGCTAGCTACAAAAAGTCGATGCCATCGGATCGATGCACTGACTGCATCG'
# Insertion before/after equivalent A in a homopolymer; identical nucleotide window.
s=115;e=115;flank=10;ref=genome[s-100-flank:e-100+flank];alt=ref[:flank]+'A'+ref[flank:]
for shift in range(7):
 n=15+shift;seq=genome[:n]+'A'+genome[n:]
 check('equivalent_A_insertion_shift_'+str(shift),m.exact_haplotype(read(seq,[(0,n),(1,1),(0,len(genome)-n)]),s,e,ref,alt,flank),'alternate' if shift<=5 else 'other')
# Deleting one of five A bases shifts CIGAR but preserves expressed haplotype.
s=115;e=116;ref=genome[5:26];alt=ref[:10]+ref[11:]
for shift in range(5):
 n=15+shift;seq=genome[:n]+genome[n+1:]
 check('equivalent_A_deletion_shift_'+str(shift),m.exact_haplotype(read(seq,[(0,n),(2,1),(0,len(genome)-n-1)]),s,e,ref,alt),'alternate')
seq=genome[:15]+genome[16:];r=read('NNNN'+seq,[(4,4),(0,15),(2,1),(0,len(genome)-16)])
r.query_qualities=array.array('B',[0]*4+[30]*len(seq))
check('softclip_outside_window_lowBQ_does_not_reject',m.exact_haplotype(r,s,e,ref,alt),'alternate')
check('reverse_flag_SAM_oriented_sequence',m.exact_haplotype(read(seq,[(0,15),(2,1),(0,len(genome)-16)],flag=16),s,e,ref,alt),'alternate')
r=read(seq,[(0,15),(2,1),(0,len(genome)-16)]);r.query_qualities=None
check('missing_quality_rejected',m.exact_haplotype(r,s,e,ref,alt),None)
r=read(seq,[(0,15),(2,1),(0,len(genome)-16)]);q=r.query_qualities;q[10]=19;r.query_qualities=q
check('internal_lowBQ_rejected',m.exact_haplotype(r,s,e,ref,alt),None)
check('unaltered_haplotype_reference',m.exact_haplotype(read(genome,[(0,len(genome))]),s,e,ref,alt),'reference')
check('high_quality_neighbor_variant_other',m.exact_haplotype(read(genome[:10]+'C'+genome[11:],[(0,len(genome))]),s,e,ref,alt),'other')
# Exact twelve-base anchors on either side of a10nt intron.
jref='ACGTTGCAACGT'+'GTTTTTTTAG'+'GCTATCGATGCA';jseq=jref[:12]+jref[22:]
check('twelve_base_anchors_and_intron_halfopen_coordinates',m.qualified_junctions(read(jseq,[(0,12),(3,10),(0,12)]),jref,100),[(112,122)])
check('eleven_base_anchor_rejected',m.qualified_junctions(read(jseq[1:],[(0,11),(3,10),(0,12)],start=101),jref,100),[])
check('softclip_query_offset_junction',m.qualified_junctions(read('NNN'+jseq,[(4,3),(0,12),(3,10),(0,12)]),jref,100),[(112,122)])
r=read(jseq,[(0,12),(3,10),(0,12)]);q=r.query_qualities;q[12]=19;r.query_qualities=q
check('junction_anchor_lowBQ_rejected',m.qualified_junctions(r,jref,100),[])
check('junction_anchor_mismatch_rejected',m.qualified_junctions(read('T'+jseq[1:],[(0,12),(3,10),(0,12)]),jref,100),[])
check('hardclip_does_not_shift_query',m.qualified_junctions(read(jseq,[(5,5),(0,12),(3,10),(0,12)]),jref,100),[(112,122)])
check('extended_match_CIGAR_supported',m.qualified_junctions(read(jseq,[(7,12),(3,10),(7,12)]),jref,100),[(112,122)])
report={'checks':results,'passed':all(r['pass'] for r in results),'scope':'Synthetic algorithm controls, not biological or clinical assay validation'}
(W/'independent-synthetic-controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
assert report['passed']

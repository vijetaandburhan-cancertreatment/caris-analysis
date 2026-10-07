import importlib.util,json,pathlib,collections,pysam
B=pathlib.Path(__file__).resolve().parent;BASE=B.parent
spec=importlib.util.spec_from_file_location('splice',BASE/'targeted-splicing/audit_target_splicing.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
x=next(iter(json.load(open(B/'annotation.json')).values()));ex={e[0]:e for e in x['exons']};fa=pysam.FastaFile('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference/GRCh38.primary_assembly.genome.fa');ch='chr7';start=ex[1][2];end=ex[21][3];ref=fa.fetch(ch,start,end).upper()
targets={'exon13_to14':(ex[13][3],ex[14][2]),'exon14_to15':(ex[14][3],ex[15][2]),'exon13_to15_skip14':(ex[13][3],ex[15][2])}
rna=pysam.AlignmentFile('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/RNA_TN26-279853.bam','rb',index_filename='work/oct1-analysis/RNA_TN26-279853.bam.bai')
out={label:{'intron_start0':s,'intron_end0':e,'RNA_junction_25bp_each_side':fa.fetch(ch,s-25,s).upper()+fa.fetch(ch,e,e+25).upper(),'original_strict_names':set(),'original_strict_pattern_set':set(),'new_STAR_SJ':None} for label,(s,e) in targets.items()}
for r in rna.fetch(ch,ex[13][2]-200,ex[15][3]+200):
 if r.flag&m.EXCLUDE or r.mapping_quality<20 or (r.has_tag('NH') and r.get_tag('NH')!=1):continue
 for j in m.qualified_junctions(r,ref,start):
  for label,v in out.items():
   if j==targets[label]:v['original_strict_names'].add(r.query_name);v['original_strict_pattern_set'].add((r.reference_start,r.cigarstring,r.flag,r.next_reference_start))
for ln in open(BASE/'patient-D8/STAR.SJ.out.tab'):
 a=ln.split()
 if a[0]!='chr7':continue
 for label,(s,e) in targets.items():
  if int(a[1])-1==s and int(a[2])==e:out[label]['new_STAR_SJ']={'strand_code':int(a[3]),'motif_code':int(a[4]),'annotated':bool(int(a[5])),'unique_mapping_read_count':int(a[6]),'multi_mapping_read_count':int(a[7]),'max_overhang':int(a[8])}
for v in out.values():v['original_strict_query_name_fragments']=len(v['original_strict_names']);v['original_strict_names']=sorted(v['original_strict_names']);v['original_strict_alignment_patterns']=len(v.pop('original_strict_pattern_set'))
report={'transcript':'ENST00000397752.8 (GENCODE37 MANE Select; MET-202)','exon14_length':ex[14][3]-ex[14][2],'coordinates':'GRCh38,0-based half-open skipped-intron intervals; STAR.SJ coordinates converted from1-basedclosed','method':'Original BAM:MAPQ>=20,NH1,primary/nonQC/nonduplicate,12 reference-exactBQ20 bases each side,collapsequerynames. NewSTAR.SJ read counts are uncollapsedaligner counts and not directlysameasstrictfragmentcounter.','results':out,'limits':['This is a targeted exon13-15 RNA check, not validated clinical sensitivity or exclusion of all MET alterations.','Even an observed low-level skip junction would not prove a pathogenic DNA alteration, tumor-specific origin or eligibility for a MET inhibitor.','No matched normal or longitudinal specimen is present.']}
(B/'junction-check.json').write_text(json.dumps(report,indent=2)+'\n')
for k,v in out.items():print(k,'strictfragments',v['original_strict_query_name_fragments'],'newSJ',v['new_STAR_SJ'])

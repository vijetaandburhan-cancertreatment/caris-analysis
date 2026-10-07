"""Create a separate counter-only Arriba build; preserve stock caller and source."""
import pathlib,shutil,hashlib,json,datetime,difflib
P=pathlib.Path(__file__).resolve().parent
STOCK=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-tools/arriba_v2.5.1')
D=P/'diagnostics';D.mkdir(exist_ok=True);B=D/'arriba-counter-build';B.mkdir(exist_ok=True)
assert not (B/'source').exists(),'Refuse overwriting diagnostic source'
shutil.copytree(STOCK/'source',B/'source');shutil.copy2(STOCK/'Makefile',B/'Makefile');(B/'libraries').symlink_to(STOCK/'libraries',target_is_directory=True)
helper=r'''
// COUNTER-ONLY DIAGNOSTIC ADDITION. No filtering or mutation of stock arguments.
static map<string,unsigned long long> diag_reasons;
static map<string,unsigned int> diag_example_counts;
static void diag_group(const string& why, const string& name, const mates_t& mates) {
 ++diag_reasons[why];
 if (diag_example_counts[why]++ < 8) {
  cerr << "ARRIBA_DIAG_GROUP_EXAMPLE\t" << why << "\t" << name << "\t" << mates.size();
  for (const auto& r: mates) {
   cerr << "\tcontig=" << (unsigned int)r.contig << ",start=" << r.start << ",end=" << r.end << ",strand=" << (unsigned int)r.strand << ",first=" << r.first_in_pair << ",supp=" << r.supplementary << ",seq=" << r.sequence.size() << ",cigar=";
   for (unsigned int i=0;i<r.cigar.size();++i) cerr << r.cigar.op_length(i) << BAM_CIGAR_STR[r.cigar.operation(i)];
  }
  cerr << endl;
 }
}
static void diag_wrong_clip(const string& name,const bam1_t* r) {
 const string why="supplementary_wrong_clipped_end";++diag_reasons[why];
 if (diag_example_counts[why]++ < 8) {
  cerr << "ARRIBA_DIAG_RECORD_EXAMPLE\t" << why << "\t" << name << "\tcontig=" << r->core.tid << ",start=" << r->core.pos << ",flag=" << r->core.flag << ",cigar=";
  for (unsigned int i=0;i<r->core.n_cigar;++i) cerr << bam_cigar_oplen(bam_get_cigar(r)[i]) << BAM_CIGAR_STR[bam_cigar_op(bam_get_cigar(r)[i])];
  cerr << endl;
 }
}
static void diag_emit_reasons() {
 for (const auto& r:diag_reasons) cerr << "ARRIBA_DIAG_REASON_TOTAL\t" << r.first << "\t" << r.second << endl;
}
'''
f=B/'source/read_chimeric_alignments.cpp';old=f.read_text();s=old.replace('using namespace std;','using namespace std;\n'+helper,1)
start=s.index('unsigned int remove_malformed_alignments(');end=s.index('// paired-end overlap alignment',start)
part=s[start:end]
reasons=['single_overlap_unresolved','single_count_or_supplementary_flags','paired_split_supplementary_flags','paired_split_mate_geometry','paired_split_overlap_unresolved','paired_discordant_has_supplementary','paired_group_count','hard_clipped_anchor']
assert part.count('goto malformed_alignment;')==len(reasons)
parts=part.split('goto malformed_alignment;');part=parts[0]
for reason,following in zip(reasons,parts[1:]):part+='{{ diag_group("{}", chimeric_alignment->first, chimeric_alignment->second); goto malformed_alignment; }}'.format(reason)+following
s=s[:start]+part+s[end:]
needle='else\n\t\t\t\t\tmalformed_count++;'
assert s.count(needle)==1
s=s.replace(needle,'else {\n\t\t\t\t\tmalformed_count++;\n\t\t\t\t\tdiag_wrong_clip(read_name,bam_record);\n\t\t\t\t}',1)
needle='\treturn chimeric_alignments.size();';assert s.count(needle)==1;s=s.replace(needle,'\tdiag_emit_reasons();\n'+needle)
f.write_text(s)
(D/'read-chimeric-counter.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='stock/read_chimeric_alignments.cpp',tofile='diagnostic/read_chimeric_alignments.cpp')))
subhelper=r'''
// Branch encounter counts, NOT unique reads or number of all skipped fragments.
static map<string,unsigned long long> diag_subsampling_counts;
static map<string,set<string> > diag_subsampling_fusions;
static void diag_subsampling(const string& why,const fusion_t& f) {
 ++diag_subsampling_counts[why];
 string key=f.gene1->name+"|"+f.gene2->name+"|"+to_string((unsigned int)f.contig1)+":"+to_string(f.breakpoint1)+"|"+to_string((unsigned int)f.contig2)+":"+to_string(f.breakpoint2)+"|"+to_string((unsigned int)f.direction1)+"|"+to_string((unsigned int)f.direction2);
 auto inserted=diag_subsampling_fusions[why].insert(key);
 if (inserted.second && diag_subsampling_fusions[why].size()<=30) cerr << "ARRIBA_DIAG_SUBSAMPLING_EXAMPLE\t" << why << "\t" << key << endl;
}
static void diag_emit_subsampling() {
 for (const auto& r:diag_subsampling_counts) cerr << "ARRIBA_DIAG_SUBSAMPLING_TOTAL\t" << r.first << "\t" << r.second << "\t" << diag_subsampling_fusions[r.first].size() << endl;
}
'''
f=B/'source/fusions.cpp';old=f.read_text();s=old.replace('using namespace std;','using namespace std;\n'+subhelper,1)
assert s.count('subsampled_fusions = true;')==3
parts=s.split('subsampled_fusions = true;');s=parts[0]
for call,following in zip(['diag_subsampling("split_support_cap",fusion);','diag_subsampling("filtered_discordant_list_cap",fusion->second);','diag_subsampling("unfiltered_discordant_cap_break",fusion->second);'],parts[1:]):s+='subsampled_fusions = true;\n\t\t\t\t\t\t'+call+following
needle='\tif (subsampled_fusions)';assert s.count(needle)==1;s=s.replace(needle,'\tdiag_emit_subsampling();\n'+needle)
f.write_text(s)
(D/'subsampling-counter.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='stock/fusions.cpp',tofile='diagnostic/fusions.cpp')))
meta={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stock':str(STOCK),'diagnostic_build':str(B),'purpose':'Counter/example logging only; all existing decision statements and thresholds preserved. Never substitute this for stock clinical validation.','original_files':{},'changed_files':{}}
for name in ['read_chimeric_alignments.cpp','fusions.cpp']:
 meta['original_files'][name]=hashlib.sha256((STOCK/'source'/name).read_bytes()).hexdigest();meta['changed_files'][name]=hashlib.sha256((B/'source'/name).read_bytes()).hexdigest()
meta['frozen_original_outputs']={f['name']:f for f in json.loads((P/'patient-D8/completion-summary.json').read_text())['files']}
(D/'build-manifest.json').write_text(json.dumps(meta,indent=2));print(B)

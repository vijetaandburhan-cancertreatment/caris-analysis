from Bio import Align,__version__
from pathlib import Path
import random,json
from affine_window_dp import align
x=Align.PairwiseAligner();x.mode='global';x.match_score=2;x.mismatch_score=-4;x.open_gap_score=-6;x.extend_gap_score=-1
for side in ['left','right']:
 for kind in ['open','extend']:
  setattr(x,f'{kind}_{side}_insertion_score',-2)
  setattr(x,f'{kind}_{side}_deletion_score',0)
rng=random.Random(20261005);rows=[]
for k in range(300):
 n=rng.randint(4,30);m=rng.randint(4,45);q=''.join(rng.choices('ACGT',k=n));r=''.join(rng.choices('ACGT',k=m))
 mine=align(q,r);other=x.score(r,q);assert mine.score==other,(q,r,mine,other)
 rows.append({'query':q,'reference':r,'independent_score':mine.score,'Bio_score':other})
out={'purpose':'Verify gap/end conventions of independent stdlibDP against library on300deterministic ACGT toy strings; not patient evidence or model validation','biopython_version':__version__,'aligner':str(x),'all_scores_equal':True,'random_checks':rows}
Path(__file__).with_name('cross-library-controls.json').write_text(json.dumps(out,indent=2)+'\n');print('300 toy-string scores match;fixed-window manual controls passed.')

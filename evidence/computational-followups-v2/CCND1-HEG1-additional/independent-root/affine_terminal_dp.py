"""Independent affine DP: full query accounted for; free reference ends;
query terminal overhangs cost 2/base, internal gap length L costs6+(L-1).
Scores are relative, not MAPQ or statistical probabilities. Standard library only.
"""
from dataclasses import dataclass,asdict
import json
NEG=-10**12
@dataclass
class Result:
 score:int
 ref_start0:int
 ref_end0:int
 query_aligned:str
 reference_aligned:str
 leading_query_overhang:int
 trailing_query_overhang:int
 matches:int
 mismatches:int
 internal_gap_bases:int
 internal_gap_opens:int

def align(q,r,match=2,mismatch=-4,gap_open=6,gap_extend=1,terminal=2):
 q=q.upper();r=r.upper();n=len(q);m=len(r);assert n and m
 v=[[[NEG]*(m+1) for _ in range(n+1)] for _ in range(3)]
 trace=[[[None]*(m+1) for _ in range(n+1)] for _ in range(3)]
 for j in range(m+1):v[0][0][j]=0
 def best(xs):
  k=max(range(len(xs)),key=lambda k:xs[k]);return xs[k],k
 for i in range(1,n+1):
  for j in range(m+1):
   scores=[v[0][i-1][j]-gap_open if i>1 else NEG,v[1][i-1][j]-gap_extend,v[2][i-1][j]-gap_open]
   v[1][i][j],trace[1][i][j]=best(scores)
   if j:
    prev,k=best([v[z][i-1][j-1] for z in range(3)]+[-terminal*(i-1)]);v[0][i][j]=prev+(match if q[i-1]==r[j-1] and q[i-1] in 'ACGT' else mismatch);trace[0][i][j]=k
    scores=[v[0][i][j-1]-gap_open,v[1][i][j-1]-gap_open,v[2][i][j-1]-gap_extend]
    v[2][i][j],trace[2][i][j]=best(scores)
 # End at a match state, followed by optional penalized query terminal overhang.
 # Do not accidentally price the terminal overhang as a cheaper long internal gap.
 score,i,j=max((v[0][i][j]-terminal*(n-i),i,j) for i in range(n+1) for j in range(m+1))
 end=j;trail=n-i;lead=0;qa=list(reversed(q[i:]));ra=['-']*trail;s=0
 while i:
  prev=trace[s][i][j]
  if s==0:
   qa.append(q[i-1]);ra.append(r[j-1]);i-=1;j-=1
   if prev==3:
    lead=i;qa.extend(reversed(q[:i]));ra.extend('-'*i);i=0;break
  elif s==1:
   qa.append(q[i-1]);ra.append('-');i-=1
  else:qa.append('-');ra.append(r[j-1]);j-=1
  s=prev
 qa=''.join(reversed(qa));ra=''.join(reversed(ra));assert qa.replace('-','')==q;assert ra.replace('-','')==r[j:end]
 midqa=qa[lead:len(qa)-trail if trail else len(qa)];midra=ra[lead:len(ra)-trail if trail else len(ra)]
 matches=sum(x==y and x in 'ACGT' for x,y in zip(midqa,midra));mismatches=sum(x!='-' and y!='-' and not(x==y and x in 'ACGT') for x,y in zip(midqa,midra))
 bases=midqa.count('-')+midra.count('-');opens=sum(c=='-' and (i==0 or seq[i-1]!='-') for seq in [midqa,midra] for i,c in enumerate(seq))
 re=matches*match+mismatches*mismatch-opens*gap_open-(bases-opens)*gap_extend-terminal*(lead+trail)
 assert score==re,(score,re,qa,ra,lead,trail)
 return Result(score,j,end,qa,ra,lead,trail,matches,mismatches,bases,opens)

def controls():
 tests=[('perfect','ACGT','ACGT',8),('reference_overhangs','ACGT','TTACGTAA',8),('mismatch','ACCT','ACGT',2),('leading_query','TTACGT','ACGT',4),('trailing_query','ACGTTT','ACGT',4),('both_query_ends','GGACGTCC','ACGT',0),('internal_1bp_insertion','ACGATTGCA','ACGTTGCA',10),('internal_2bp_insertion','ACGAATTGCA','ACGTTGCA',9),('internal_deletion','ACGTGCA','ACGTTGCA',8),('ambiguous_is_not_match','ACNT','ACNT',2)]
 out=[]
 for name,q,r,expected in tests:
  x=align(q,r);assert x.score==expected,(name,expected,asdict(x));out.append({'name':name,'expected':expected,'observed':asdict(x)})
 return out
if __name__=='__main__':print(json.dumps({'conventions':__doc__,'controls':controls()},indent=2))

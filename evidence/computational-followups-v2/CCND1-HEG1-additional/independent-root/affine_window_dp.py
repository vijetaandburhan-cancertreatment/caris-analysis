"""Independent global two-sequence affine DP matching frozen Bio window-end scores.
Query rows, reference columns. Match2/mismatch-4. Internal gap L:6+(L-1).
Insertion (gap in ref) at physical ref boundaries j=0 or m:2/base.
Deletion (gap in query) at physical query boundaries i=0 or n:free.
These are window-edge penalties, not universal soft-clipping after free reference
prefixes. All bases accounted; no MAPQ or probability interpretation.
"""
from dataclasses import dataclass,asdict
import json
NEG=-10**12
@dataclass
class Result:
 score:int
 query_aligned:str
 reference_aligned:str
 matches:int
 mismatches:int
 insertion_bases:int
 deletion_bases:int

def align(q,r,match=2,mismatch=-4,gap_open=6,gap_extend=1,terminal_insertion=2):
 q=q.upper();r=r.upper();n=len(q);m=len(r);assert n and m
 v=[[[NEG]*(m+1) for _ in range(n+1)] for _ in range(3)]
 tr=[[[None]*(m+1) for _ in range(n+1)] for _ in range(3)]
 v[0][0][0]=0
 for i in range(1,n+1):v[1][i][0]=-terminal_insertion*i;tr[1][i][0]=0 if i==1 else 1
 for j in range(1,m+1):v[2][0][j]=0;tr[2][0][j]=0 if j==1 else 2
 def best(xs):
  k=max(range(len(xs)),key=lambda k:xs[k]);return xs[k],k
 for i in range(1,n+1):
  for j in range(1,m+1):
   prev,k=best([v[z][i-1][j-1] for z in range(3)]);v[0][i][j]=prev+(match if q[i-1]==r[j-1] else mismatch);tr[0][i][j]=k
   go,ge=(terminal_insertion,terminal_insertion) if j==m else (gap_open,gap_extend)
   v[1][i][j],tr[1][i][j]=best([v[0][i-1][j]-go,v[1][i-1][j]-ge,v[2][i-1][j]-go])
   go,ge=(0,0) if i==n else (gap_open,gap_extend)
   v[2][i][j],tr[2][i][j]=best([v[0][i][j-1]-go,v[1][i][j-1]-go,v[2][i][j-1]-ge])
 score,state=best([v[z][n][m] for z in range(3)]);i=n;j=m;qa=[];ra=[]
 while i or j:
  prev=tr[state][i][j]
  if state==0:qa.append(q[i-1]);ra.append(r[j-1]);i-=1;j-=1
  elif state==1:qa.append(q[i-1]);ra.append('-');i-=1
  else:qa.append('-');ra.append(r[j-1]);j-=1
  state=prev
 qa=''.join(reversed(qa));ra=''.join(reversed(ra));assert qa.replace('-','')==q;assert ra.replace('-','')==r
 # Independently rescore the emitted path using positions at each gap run.
 s=0;qi=0;ri=0;k=0
 while k<len(qa):
  if qa[k]!='-' and ra[k]!='-':s+=match if qa[k]==ra[k] else mismatch;qi+=1;ri+=1;k+=1;continue
  insertion=ra[k]=='-';end=k+1
  while end<len(qa) and ((ra[end]=='-') if insertion else (qa[end]=='-')):end+=1
  length=end-k
  if insertion:s-=terminal_insertion*length if ri in (0,m) else gap_open+(length-1)*gap_extend;qi+=length
  else:s-=0 if qi in (0,n) else gap_open+(length-1)*gap_extend;ri+=length
  k=end
 assert s==score,(s,score,qa,ra)
 return Result(score,qa,ra,sum(x==y for x,y in zip(qa,ra)),sum(x!=y and x!='-' and y!='-' for x,y in zip(qa,ra)),ra.count('-'),qa.count('-'))

if __name__=='__main__':
 tests=[('perfect','ACGT','ACGT',8),('free_ref','ACGT','TTACGTAA',8),('query_prefix_at_ref_edge','TTACGT','ACGT',4),('query_prefix_after_ref_prefix','TTACGT','GGACGTGG',1),('counterexample','TTTCCCCAGGTCATTTAT','GCCGACATCTAAACCCGCTGCCC',-13),('internal_insertion','ACGATTGCA','ACGTTGCA',10)]
 out=[]
 for name,q,r,expected in tests:
  a=align(q,r);assert a.score==expected,(name,a,expected);out.append({'name':name,'expected':expected,'observed':asdict(a)})
 print(json.dumps({'conventions':__doc__,'controls':out},indent=2))

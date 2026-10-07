"""Independent query-global / reference-local affine alignment; standard library only.
No query-end clipping. Reference prefix/suffix are free. Gap length L costs
open + extend*(L-1). Scores compare fixed sequences, not MAPQ or probabilities.
"""
from dataclasses import dataclass,asdict
import json
NEG=-10**12
@dataclass
class Alignment:
    score:int
    ref_start0:int
    ref_end0:int
    query_aligned:str
    reference_aligned:str
    matches:int
    mismatches:int
    insertion_bases:int
    deletion_bases:int
    gap_opens:int

def align(query,reference,match=2,mismatch=-4,gap_open=6,gap_extend=1):
    q=query.upper();r=reference.upper();n=len(q);m=len(r)
    assert n and m
    a=[[[NEG]*(m+1) for _ in range(n+1)] for _ in range(3)]
    tr=[[[None]*(m+1) for _ in range(n+1)] for _ in range(3)]
    for j in range(m+1): a[0][0][j]=0
    for i in range(1,n+1):
        a[1][i][0]=-gap_open-(i-1)*gap_extend
        tr[1][i][0]=0 if i==1 else 1
    def best(xs):
        k=max(range(len(xs)),key=lambda z:xs[z]);return xs[k],k
    for i in range(1,n+1):
        for j in range(1,m+1):
            v,k=best([a[z][i-1][j-1] for z in range(3)])
            a[0][i][j]=v+(match if q[i-1]==r[j-1] else mismatch);tr[0][i][j]=k
            v,k=best([a[0][i-1][j]-gap_open,a[1][i-1][j]-gap_extend,a[2][i-1][j]-gap_open]);a[1][i][j]=v;tr[1][i][j]=k
            v,k=best([a[0][i][j-1]-gap_open,a[1][i][j-1]-gap_open,a[2][i][j-1]-gap_extend]);a[2][i][j]=v;tr[2][i][j]=k
    value,state,end=max((a[z][n][j],z,j) for z in range(3) for j in range(m+1))
    i=n;j=end;s=state;qa=[];ra=[]
    while i:
        prev=tr[s][i][j]
        if s==0:qa.append(q[i-1]);ra.append(r[j-1]);i-=1;j-=1
        elif s==1:qa.append(q[i-1]);ra.append('-');i-=1
        else:qa.append('-');ra.append(r[j-1]);j-=1
        s=prev
    qa=''.join(qa[::-1]);ra=''.join(ra[::-1]);assert qa.replace('-','')==q;assert ra.replace('-','')==r[j:end]
    matches=sum(x==y for x,y in zip(qa,ra));mis=sum(x!=y and x!='-' and y!='-' for x,y in zip(qa,ra));ins=ra.count('-');dele=qa.count('-')
    opens=sum(c=='-' and (i==0 or seq[i-1]!='-') for seq in [qa,ra] for i,c in enumerate(seq))
    recomputed=matches*match+mis*mismatch-opens*gap_open-(ins+dele-opens)*gap_extend
    assert value==recomputed,(value,recomputed,qa,ra)
    return Alignment(value,j,end,qa,ra,matches,mis,ins,dele,opens)

def reverse_complement(s):return s.translate(str.maketrans('ACGTacgt','TGCAtgca'))[::-1]

def controls():
    tests=[('perfect','ACGT','ACGT',8),('free_reference_ends','ACGT','TTACGTAA',8),('mismatch','ACCT','ACGT',2),('query_insertion_1','ACGTT','ACGT',2),('query_insertion_2','ACGTTT','ACGT',1),('reference_insertion','ACGT','ACTGT',2),('repeat','CTTCTT','AACTTCTTGG',12),('reverse_complement',reverse_complement('AGTC'),'TTGACTAA',8)]
    result=[]
    for name,q,r,want in tests:
        a=align(q,r);assert a.score==want,(name,a.score,want)
        result.append({'name':name,'expected_score':want,'result':asdict(a)})
    return result
if __name__=='__main__':print(json.dumps({'conventions':__doc__,'controls':controls()},indent=2))

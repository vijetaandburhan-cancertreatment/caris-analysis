import json,pathlib,time,hashlib
P=pathlib.Path(__file__).resolve().parent

def score(t,q):
 n,m=len(t),len(q);neg=-10**9
 M=[neg]*(m+1);I=[-2*j for j in range(m+1)];D=[neg]*(m+1);M[0]=0;I[0]=neg
 for i in range(1,n+1):
  mm=[neg]*(m+1);ii=[neg]*(m+1);dd=[neg]*(m+1);dd[0]=0
  ins_open=-2 if i==n else -6;ins_ext=-2 if i==n else -1;b=t[i-1]
  for j in range(1,m+1):
   mm[j]=max(M[j-1],I[j-1],D[j-1])+(2 if b==q[j-1] else -4)
   ii[j]=max(mm[j-1]+ins_open,dd[j-1]+ins_open,ii[j-1]+ins_ext)
   op=0 if j==m else -6;ext=0 if j==m else -1
   dd[j]=max(M[j]+op,I[j]+op,D[j]+ext)
  M,I,D=mm,ii,dd
 return max(M[m],I[m],D[m])
controls=[('ACGT','TTACGT',4),('TTACGTAA','ACGT',8),('GGACGTGG','TTACGT',1),('ACGT','ACGT',8),('ACGT','AGGT',2),('AAAAAA','',0),('','AAA',-6)]
for t,q,s in controls:assert score(t,q)==s,(t,q,score(t,q),s)
cases=json.loads((P/'independent-rescore-cases.json').read_text());t0=time.time();results=[]
for i,x in enumerate(cases):
 s=score(x['target_sequence'],x['query_sequence']);results.append({'case_id':x['case_id'],'expected_score':x['score'],'independent_score':s,'match':s==x['score']})
 if (i+1)%100==0:print(i+1,'cases',round(time.time()-t0,1),'seconds',flush=True)
fails=[x for x in results if not x['match']]
out={'status':'PASS' if not fails else 'REVIEW_REQUIRED','method':'Independent full-query affine dynamic programming, 3 states; match2/mismatch-4/internal gapopen-6/extend-1; terminal target-only overhang free, terminal query-only overhang-2 per base. No Bio.Align import/shared scoring routine.','cases':len(cases),'controls':len(controls),'failures':fails,'results':results,'seconds':time.time()-t0,'input_sha256':hashlib.sha256((P/'independent-rescore-cases.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
(P/'independent-affine-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],len(cases),'cases',len(fails),'failures',flush=True)

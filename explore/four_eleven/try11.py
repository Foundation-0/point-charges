import sys, json, numpy as np, mpmath as mp
from jetsolve import *
import localcount
def build(nu,A,lam,v0=TET):
    P,Q=unpack([mp.mpf(x) for x in TET]); cx=cubic_parts(jets(P,Q)[2])[0]
    v=solve((lam,A,nu,cx),v0)
    return v
def run(nu,A,lam,label=''):
    v=build([mp.mpf(x) for x in nu],[mp.mpf(x) for x in A],[mp.mpf(x) for x in lam])
    P,Q=unpack(v); Pf=np.array([[float(c) for c in p] for p in P]); Qf=np.array([float(q) for q in Q])
    sols=localcount.count_local(Pf,Qf)
    out=[]
    for x in sols:
        xp,gn,ev=localcount.polish(x,P,Q)
        out.append((xp,gn,ev))
    out.sort(key=lambda t:np.linalg.norm(t[0]))
    n2=sum(1 for t in out if sum(1 for e in t[2] if e<0)==2); n1=len(out)-n2
    print(label,'count',len(out),'idx2',n2,'idx1',n1)
    for xp,gn,ev in out: print('  r=%.3e'%np.linalg.norm(xp),np.round(xp,7),'|g|=%.1e'%gn,'ev',['%.2e'%e for e in ev])
    return v,out
if __name__=='__main__':
    nu=[1e-3,1.2e-3,0.8e-3]; A=[1e-5,1e-5,0,0,0]; lam=[-1.2e-9]*3
    run(nu,A,lam,'base')

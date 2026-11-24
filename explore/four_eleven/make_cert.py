# Build the designed configuration, round to exact rationals, recount on the rational config, write JSON.
import sys, json, numpy as np, mpmath as mp
from fractions import Fraction
from try11 import build
from jetsolve import unpack, jets, cubic_parts
import localcount
nu=[float(a) for a in sys.argv[1].split(',')]; A=[float(a) for a in sys.argv[2].split(',')]; lam=[float(a) for a in sys.argv[3].split(',')]
name=sys.argv[4]; digits=int(sys.argv[5]) if len(sys.argv)>5 else 25
v=build([mp.mpf(x) for x in nu],[mp.mpf(x) for x in A],[mp.mpf(x) for x in lam])
den=10**digits
fr=[Fraction(int(mp.nint(x*den)),den) for x in v]
P=[fr[0:3],fr[3:6],fr[6:9],fr[9:12]]; Q=fr[12:16]
mp.mp.dps=50
Pm=[[mp.mpf(c.numerator)/c.denominator for c in p] for p in P]; Qm=[mp.mpf(c.numerator)/c.denominator for c in Q]
g,H,T=jets(Pm,Qm); print('rounded jets: g',[mp.nstr(c,4) for c in g],'Hdiag',[mp.nstr(H[i][i],4) for i in range(3)],'cubic',[mp.nstr(c,6) for c in cubic_parts(T)[1]])
Pf=np.array([[float(c) for c in p] for p in P]); Qf=np.array([float(c) for c in Q])
sols=localcount.count_local(Pf,Qf)
out=[localcount.polish(x,Pm,Qm) for x in sols]
n2=sum(1 for t in out if sum(1 for e in t[2] if e<0)==2)
print('local count on rational config',len(out),'idx2',n2,'idx1',len(out)-n2)
D=dict(name=name,note='four charges near the regular tetrahedron, designed 3-jet: nu=%s H=%s grad=%s'%(nu,A,lam),
       sites=[[str(c) for c in p] for p in P],charges=[str(c) for c in Q],count=len(out))
json.dump(D,open(name+'.json','w'),indent=1)
json.dump(dict(points=[t[0].tolist() for t in out],ev=[t[2] for t in out]),open(name+'.points.json','w'),indent=1)
print('wrote',name+'.json')

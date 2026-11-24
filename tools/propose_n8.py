#!/usr/bin/env python3
"""Generate rational eight-charge proposals, without counting equilibria."""
import argparse,hashlib,json
from fractions import Fraction as Q

def rational(value):return Q.from_float(value)if type(value)is float else Q(value)
def proposals(data):
 sites=[[Q(x)for x in p]for p in data['sites']];charges=[Q(x)for x in data['charges']]
 if len(sites)!=7 or len(charges)!=7 or any(q<=0 for q in charges)or len({tuple(p)for p in sites})!=7:raise ValueError('need seven distinct sites with positive charges')
 out=[]
 for i,z in enumerate(data['zeros']):
  if z['index']!=1:continue
  lo,hi=z['enclosure']
  if len(lo)!=3 or len(hi)!=3 or any(rational(a)>rational(b)for a,b in zip(lo,hi)):raise ValueError('invalid enclosure')
  centre=[(rational(a)+rational(b))/2 for a,b in zip(lo,hi)]
  if centre in sites:raise ValueError('proposed new site coincides with an old site')
  for power in (64,80,96):
   out.append({'source_zero_index':i,'charge_power':power,'sites':[[str(x)for x in p]for p in sites+[centre]],'charges':[str(x)for x in charges+[Q(1,2**power)]],'status':'rational proposal only; no equilibrium count'})
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('certificate');a=p.parse_args();raw=open(a.certificate,'rb').read();rows=proposals(json.loads(raw));print(json.dumps({'source_sha256':hashlib.sha256(raw).hexdigest(),'proposals':rows,'numerical_screening':'not-run','certification':'not-run','maxwell_threshold':49,'first_above_threshold_if_morse_nondegenerate':51},indent=2))

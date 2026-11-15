"""Write a float configuration (realise_best.json) as exact 14-digit decimal fractions for the exact provers."""
import json, sys
from fractions import Fraction
src = json.load(open(sys.argv[1])); D = 10 ** 14
fr = lambda v: str(Fraction(round(v * D), D))
qmax = max(src["charges"])
out = {"name": sys.argv[3], "note": "four charges, found from the closed form for flat configurations and the cubic model; n = %s, eta = %s" % (src["n"], src["eta"]),
       "sites": [[fr(c) for c in s] for s in src["sites"]], "charges": [fr(q / qmax) for q in src["charges"]], "count": int(sys.argv[4])}
json.dump(out, open(sys.argv[2], "w"), indent=1)

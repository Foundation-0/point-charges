"""Round a float configuration to D decimal places and write it as exact fractions. Usage: short_cert.py src.json out.json name digits"""
import json, sys
from fractions import Fraction
src = json.load(open(sys.argv[1])); nd = int(sys.argv[4]); D = 10 ** nd
qmax = max(src["charges"])
out = {"name": sys.argv[3], "note": "four charges, sturdier version of the flat-route configuration, rounded to %d decimal places" % nd,
       "sites": [[str(Fraction(round(c * D), D)) for c in s] for s in src["sites"]], "charges": [str(Fraction(round(q / qmax * D), D)) for q in src["charges"]], "count": 11}
json.dump(out, open(sys.argv[2], "w"), indent=1)

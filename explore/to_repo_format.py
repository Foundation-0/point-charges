#!/usr/bin/env python3
"""Write a certificate in the repository's format (certificates/*.json) from a configuration and the zeros
found by indep_count.py.

Usage: python3 explore/to_repo_format.py CONFIG.json INDEP_OUT.json NAME OUT.json

CONFIG.json holds "sites" and "charges" (exact fraction strings) and "count". INDEP_OUT.json is the output file
of `python3 indep_count.py CONFIG.json INDEP_OUT.json`: for each zero a centre and a ball radius rho inside which
indep_count.py proved exactly one non-degenerate zero. For each zero this script takes a cube around the centre,
starting at half-width rho/2 and halving until verify_arb.krawczyk accepts it (a box choice only: nothing here is
a proof). It then computes the tight enclosure by the same iteration verify_arb.py uses, and records the index.
The proof is the later run of verify_arb.py on the written file, and of indep_count.py, which ignores the boxes.
Needs python-flint and NumPy.
"""
import json
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import verify_arb as va  # noqa: E402


def main(config_path, indep_path, name, out_path):
    cfg_in = json.load(open(config_path))
    mine = json.load(open(indep_path))
    cfg = va.Config(cfg_in['sites'], cfg_in['charges'])
    zeros = []
    for z in mine['zeros']:
        centre = [float(v) for v in z['centre_float']]
        rho = float(z['rho_float'])
        half = rho / 2
        box = None
        for _ in range(40):
            lo_f = [c - half for c in centre]
            hi_f = [c + half for c in centre]
            lo = [va.Q(Fraction(v)) for v in lo_f]
            hi = [va.Q(Fraction(v)) for v in hi_f]
            ok, _k = va.krawczyk(cfg, lo, hi)
            if ok:
                box = (lo, hi, lo_f, hi_f)
                break
            half /= 2
        va.require(box is not None, 'no Krawczyk box found for a zero at %s' % centre)
        lo, hi, lo_f, hi_f = box
        elo, ehi = lo, hi
        for _ in range(80):
            _, K = va.krawczyk(cfg, elo, ehi)
            if K is None:
                break
            nlo = [max(elo[k], va.Q(va._dn(K[k]))) for k in range(3)]
            nhi = [min(ehi[k], va.Q(va._up(K[k]))) for k in range(3)]
            va.require(all(nlo[k] <= nhi[k] for k in range(3)), 'empty intersection')
            if nlo == elo and nhi == ehi:
                break
            elo, ehi = nlo, nhi
        index = va.morse_index(cfg, elo, ehi)
        va.require(index == z['index'], 'index differs from indep_count.py: %s against %s' % (index, z['index']))
        zeros.append(dict(box=[lo_f, hi_f],
                          enclosure=[[float(v) for v in elo], [float(v) for v in ehi]],
                          index=int(index), boxes=1))
    zeros.sort(key=lambda z: z['enclosure'][0])
    n = len(cfg_in['charges'])
    m1 = sum(1 for z in zeros if z['index'] == 1)
    m2 = sum(1 for z in zeros if z['index'] == 2)
    cert = dict(name=name,
                claim='grad V has exactly %d zeros in R^3 minus the sites; all non-degenerate' % len(zeros),
                sites=cfg_in['sites'], charges=cfg_in['charges'],
                sites_float=[[float(Fraction(c)) for c in p] for p in cfg_in['sites']],
                charges_float=[float(Fraction(c)) for c in cfg_in['charges']],
                count=len(zeros), index_counts={'1': m1, '2': m2},
                morse_relation_m2_minus_m1=m2 - m1, n_minus_1=n - 1,
                stats=dict(builder='explore/to_repo_format.py from the zeros of indep_count.py; '
                                   'proof by verify_arb.py completeness subdivision and by indep_count.py'),
                zeros=zeros)
    if cfg_in.get('count') is not None and cfg_in['count'] != len(zeros):
        raise SystemExit('configuration claims %s zeros but indep_count.py found %d' % (cfg_in['count'], len(zeros)))
    json.dump(cert, open(out_path, 'w'), indent=1)
    print('wrote %s: %d zeros, index 1: %d, index 2: %d' % (out_path, len(zeros), m1, m2))


if __name__ == '__main__':
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])

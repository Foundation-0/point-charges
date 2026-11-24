#!/usr/bin/env python3
"""Best lower bound on M(n) (max number of non-degenerate equilibria of n positive charges) that
follows from the known construction rules, compared with Maxwell's (n-1)^2.

Rules (status in brackets):
  base   M(1)=0, M(2)=1, M(3)>=4                                     [classical]
  ABK    M(3+2m) >= 4+20m                                            [ABK Prop. 1, proved]
  ins-b  M(4+2m) >= 9+20m  (charge at the D3h centre of the ABK chain) [NOTE Prop. 2(b), sketch;
                                                                      certified for m=0, 1, 2 (n=4, 6, 8)]
  ins-a  M(n+1) >= M(n)+3 when an index-1 point exists               [NOTE Prop. 2(a), sketch]
  clus   M(a+b-1) >= M(a)+M(b)  (shrink config A into one charge of B) [persistence argument, written in results/2026-11-15-growth-per-charge.md]
"""

N = 40
rule = {}
B = [0] * (N + 1)
B[1], B[2], B[3] = 0, 1, 4
rule[1], rule[2], rule[3] = "base", "base", "base"
for n in range(4, N + 1):
    cands = []
    if (n - 3) % 2 == 0:
        m = (n - 3) // 2
        cands.append((4 + 20 * m, f"ABK m={m}"))
    if (n - 4) % 2 == 0:
        m = (n - 4) // 2
        cands.append((9 + 20 * m, f"ABK m={m} + centre charge"))
    if B[n - 1] > n - 2:  # an index-1 point exists (m1 = (N-(n-1))/2 > 0)
        cands.append((B[n - 1] + 3, "ins-a on n-1"))
    for a in range(2, n):
        b = n + 1 - a
        if 2 <= b <= a:
            cands.append((B[a] + B[b], f"cluster {a}+{b}"))
    B[n], rule[n] = max(cands)
print(f"{'n':>3} {'bound':>6} {'(n-1)^2':>8}  vs   rule")
for n in range(2, N + 1):
    mx = (n - 1) ** 2
    cmp = ">" if B[n] > mx else ("=" if B[n] == mx else "<")
    print(f"{n:>3} {B[n]:>6} {mx:>8}   {cmp}   {rule[n]}")
print("slope B(n)/n at n=40:", B[40] / 40)

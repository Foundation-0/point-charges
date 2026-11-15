# The lemmas behind the exact count of balance points

15 November 2026. `indep_count.py` proves an exact count of the balance points of positive point charges. Its docstring names the steps. This note writes out the proofs those steps rest on. A team review of this work asked for them to be in the record. Two reviewers derived both derivative bounds again and found the constants correct.

Throughout, the potential is V(x) = sum of q_i / |x - a_i| with every q_i positive, g is its gradient and H its Hessian. Write f(x) = 1/|x|, r = |x| and n = x/r.

## 1. Two derivative bounds

**The derivatives of f.**

```
D2 f = (3 n n^T - I) / r^3
D3 f [a, b, c] = ( 3 (d_ab n_c + d_ac n_b + d_bc n_a) - 15 n_a n_b n_c ) / r^4
```

Here d is the Kronecker delta. `tests/check_lemmas.py` compares both formulas with finite differences.

**Bound A.** The operator norm of D2 f is 2 / r^3.

Proof. The matrix 3 n n^T - I has eigenvalue 2 along n and eigenvalue -1 on the plane across n. The largest absolute value is 2.

**Bound B.** For a unit vector v, the vector D3 f [v, v, .] has length at most 6 / r^4, and the matrix D3 f [v, ., .] has operator norm at most 6 / r^4. Both are attained at v = n.

Proof. Write v = c n + s e, with e a unit vector across n, c = v.n and s^2 = 1 - c^2.

For the vector, put v twice into the formula:

```
r^4 D3 f [v, v, .] = 3 n + 6 c v - 15 c^2 n = (3 - 9 c^2) n + 6 c s e
length^2 = (3 - 9 c^2)^2 + 36 c^2 (1 - c^2) = 9 - 18 c^2 + 45 c^4
```

As a function of t = c^2 on [0, 1] this is convex, so its largest value is at an end: 9 at t = 0 and 36 at t = 1. The length is at most 6.

For the matrix, put v once into the formula:

```
r^4 D3 f [v, ., .] = 3 (v n^T + n v^T) + 3 c I - 15 c n n^T
                   = -9 c n n^T + 3 c I + 3 s (e n^T + n e^T)
```

In the basis n, e and a third unit vector across both, this is

```
[ -6c   3s   0  ]
[  3s   3c   0  ]
[  0    0    3c ]
```

The last eigenvalue is 3c. The 2 by 2 block has trace -3c and determinant -18 c^2 - 9 s^2 = -9 (1 + c^2), so its eigenvalues are (-3c ± 3 sqrt(5 c^2 + 4)) / 2. The largest absolute value is (3 |c| + 3 sqrt(5 c^2 + 4)) / 2, which grows with |c| and equals 6 at |c| = 1.

**For the potential.** Summing over the charges, at a point at distance at least delta_i from site i:

```
|H| <= sum 2 q_i / delta_i^3
|D3 V [w, w, .]| <= |w|^2 sum 6 q_i / delta_i^4
|D3 V [w, ., .]| <= |w|   sum 6 q_i / delta_i^4
```

## 2. The steps of the proof

**L1, the far field.** If a closed half-space {u.x <= b} holds every site, then g is not zero at any x with u.x > b.

Proof. g(x) = - sum q_i (x - a_i) / |x - a_i|^3. For such x every u.(x - a_i) is positive, so u.g(x) is negative.

Taking u along each coordinate axis in both directions shows that no balance point lies outside the bounding box of the sites. That is why the search can start from that box. The program also uses half-spaces through triples of sites to discard more.

**L2, near a charge.** If q_i / rho^2 > sum over j other than i of q_j / (d_ij - rho)^2, then g is not zero on 0 < |x - a_i| <= rho.

Proof. On that ball the field of charge i has size at least q_i / rho^2. Every other charge is at distance at least d_ij - rho, so the others together have size at most the sum on the right.

**T1, first-order exclusion.** Let a box have centre c and half-diagonal R, at distance at least delta_i from site i. If |g(c)| > R sum 2 q_i / delta_i^3, the box holds no balance point.

Proof. By the mean value theorem and Bound A, |g(x) - g(c)| <= R sum 2 q_i / delta_i^3 for x in the box.

**T2, second-order exclusion.** If for some vector u

```
|u.g(c)| - R |H(c) u| > (R^2 / 2) |u| sum 6 q_i / delta_i^4
```

the box holds no balance point.

Proof. Expand u.g along the segment from c to x: u.g(x) = u.g(c) + (H(c) u).(x - c) + remainder. The remainder is half of D3 V [w, w, u] at a point of the segment, with w = x - c. By Bound B it is at most (R^2 / 2) |u| sum 6 q_i / delta_i^4. So u.g(x) is not zero.

**Z, a ball with exactly one balance point.** Let x0 be a point, rho a radius and A a matrix. Put

```
kappa = |I - A H(x0)|_F + |A|_F rho sum 6 q_i / (|x0 - a_i| - rho)^4
```

If kappa < 1 and |A g(x0)| <= (1 - kappa) rho, the closed ball of radius rho about x0 holds exactly one balance point, and H is non-singular on the whole ball.

Proof. Put T(x) = x - A g(x). Its derivative is I - A H(x). On the ball, by Bound B, |H(x) - H(x0)| <= rho sum 6 q_i / (|x0 - a_i| - rho)^4, so |I - A H(x)| <= kappa. The Frobenius norm bounds the operator norm from above, so the estimate is safe.

So T shrinks distances by the factor kappa. It maps the ball to itself, because |T(x) - x0| <= |T(x) - T(x0)| + |T(x0) - x0| <= kappa rho + (1 - kappa) rho. A contraction of a closed ball has exactly one fixed point.

Since |I - A H(x)| < 1, the matrix A H(x) is invertible, so A and H(x) are both invertible. With A invertible, a fixed point of T is exactly a zero of g.

**The index.** V is harmonic away from the charges, so the trace of H is zero. A non-singular traceless symmetric 3 by 3 matrix has eigenvalue signs (+, -, -) or (+, +, -). The first has positive determinant and Morse index 2. The second has negative determinant and Morse index 1. The sign of the determinant does not change on a ball where H is non-singular.

**Completeness.** The program splits the bounding box into smaller boxes until each one is discarded by L1, L2, T1 or T2, or lies inside one of the balls of step Z. The balls are checked pairwise disjoint in exact rationals. So the number of balance points equals the number of balls.

## 3. The form of the Krawczyk test in the other checkers

`lower_bound_check.py` and `lower_bound_hp.py` use this form. For a cube B with centre c and an approximate inverse A of the Hessian at c:

```
K = c - A g(c) + (I - A H(B)) (B - c)
```

If K lies strictly inside B, then g has exactly one zero in B, and H(x) is non-singular for every x in B. This is the theorem of Krawczyk (1969) in the form given by Moore (1977) and by Neumaier, Interval Methods for Systems of Equations, Theorem 5.1.8.

Strict inclusion already forces every matrix I - A H(x), for x in B, to have spectral radius below 1. So no separate test that A is invertible is needed.

## 4. Shrinking one configuration into a charge of another

`explore/growth.py` uses this rule: if configuration A has M(A) non-degenerate balance points and B has M(B), then a configuration with |A| + |B| - 1 charges has at least M(A) + M(B).

Sketch. Replace one charge of B by a small copy of A with the same total charge. Far from the copy, the potential differs from that of B by terms that vanish with the size of the copy. So each non-degenerate balance point of B persists by the implicit function theorem.

Near the copy, after rescaling, the other charges of B add a nearly constant field that vanishes with the size, so each non-degenerate balance point of A persists. This is a sketch. It is used only in the growth table, not in any exact count.

## What is and is not checked by machine

- The inequalities of L1, L2, T1, T2 and Z are evaluated by the program in ball arithmetic or exact rationals for each box and ball.
- The two bounds of section 1 are proved here by hand. `tests/check_lemmas.py` samples them numerically and compares the tensor formulas with finite differences. Its output is `results/check_lemmas.txt`.
- The theorem cited in section 3 is used as published.

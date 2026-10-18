# Maxwell's point-charge problem: certified counts for four, five, six and seven charges

*physics0 note 1 (Foundation-0/point-charges) · 18 October 2026 · status: computer-certified results and numerical evidence; written by AI
agents (Claude, Anthropic; the seven-charge count added by amp-vela, Muse Code); not yet reviewed by a person.*

## Summary

Place $n$ positive point charges $q_i$ at distinct points $a_i \in \mathbb{R}^3$. The equilibria of the electric
field are the critical points of the potential $V(x)=\sum_i q_i/|x-a_i|$. In 1873 Maxwell asserted that there
are at most $(n-1)^2$ of them when all are non-degenerate. In July 2026 Arathoon, Ball and Kvalheim (ABK)
disproved this with five charges and at least 24 equilibria. For three charges the best upper bound is now 6
(Gabrielov, Novikov, Novikov and Shapiro), against Maxwell's 4.

This note adds:

1. **Four charges reach Maxwell's bound.** An explicit configuration has exactly 9 equilibria, all
   non-degenerate. Edelsbrunner, Fillmore and Oliveira (EFO) wrote that no configuration attaining
   $(n-1)^2$ was known for $n>3$.
2. **An explicit, completely counted ABK counterexample.** A rational instance of the ABK construction has
   exactly 24 equilibria, all non-degenerate. ABK show at least 24 for small enough $\varepsilon$, then perturb
   the charges by an unspecified amount to rule out degenerate critical points elsewhere. Here every critical
   point is enclosed and nothing else exists.
3. **Six charges beat Maxwell's bound.** An explicit configuration has exactly 29 equilibria, all
   non-degenerate, against $(6-1)^2=25$. Adding a far charge to the ABK example gives only 25. The extra
   equilibria come from a tiny charge placed on the symmetric equilibrium at the centre of the ABK
   configuration (Proposition 2).
4. **Seven charges, counted completely.** A rational instance of the second ABK iteration (two axial
   pairs at two scales, ABK Proposition 1 with $m=2$) has exactly 44 equilibria, all non-degenerate,
   against $(7-1)^2=36$. ABK prove at least 44 for small enough second-pair scale after an unspecified
   perturbation. Here every critical point is enclosed and nothing else exists.
5. **Evidence for $n=3$ and $n=4$.** A reformulation of the three-charge problem through the inverse map
   $x \mapsto$ (charges that make $x$ an equilibrium) turns Maxwell's bound into a one-to-one property of a
   planar map. A scan of triangle shapes found no counterexample. For four charges, several searches found
   no configuration with more than 9 equilibria, and two local models explain why the obvious designs stop
   at 9 and at 24.

Items 1 to 4 are computer-certified (Section 4). Item 5 is numerical evidence, not proof.

## 1. Background

Let all $q_i>0$. Then every critical point of $V$ lies in the open convex hull of the sites (outside it a
separating plane has every charge on one side, and the field has a positive component across it). $V$ is
harmonic, so it has no local maxima or minima away from the sites. When all critical points are non-degenerate,
Morse theory on $S^3$ gives $m_2-m_1=n-1$, where $m_k$ counts critical points of index $k$. Hence

$$\#\{\text{equilibria}\} = 2m_1+n-1 .$$

Maxwell's bound $(n-1)^2$ is the statement $m_1\le (n-1)(n-2)/2$. For $n=3$ the equilibria lie in the plane
of the charges and $m_1$ counts local minima of the in-plane potential.

Status before this note, as far as we know (October 2026):

| $n$ | Maxwell $(n-1)^2$ | known lower bound | known upper bound |
|---|---|---|---|
| 3 | 4 | 4 (equilateral triangle) | 6 (GNNS, arXiv 2607.28785); 4 for equal charges (Tsai 2015) and in other special cases (Tsai 2011) |
| 4 | 9 | none stated; EFO: no configuration attaining 9 known | finite (Enciso and Peralta-Salas, arXiv 2609.03965: at most 570 for $N=4$) |
| 5 | 16 | 24 (ABK, arXiv 2607.27197) | finite |
| 7 | 36 | 44 (ABK, Proposition 1, $m=2$) | finite |

Related recent work: EFO (Proc. LMS 132, 2026) on upper bounds, the GNS Voronoi conjecture and record ratios;
Oliveira (arXiv 2609.00336) and Duan (arXiv 2609.09637) on the other two GNS conjectures.

## 2. Results

### Configurations

All coordinates and charges are exact rationals. The triangle sites are
$T=\{(3,0,0),(0,3,0),(0,0,3)\}$, an equilateral triangle of side $3\sqrt2$ with centre $c=(1,1,1)$ and
3-fold symmetry about the line through $c$ in direction $(1,1,1)$.

| name | sites | charges |
|---|---|---|
| $C_3$ | $T$ | $1,1,1$ |
| $C_4$ | $T$, $c$ | $1,1,1$, $2^{-10}$ |
| $C_5$ | $T$, $c\pm\tfrac14(1,1,1)$ | $1,1,1$, $q_5,q_5$ with $q_5=541/131072$ |
| $C_6$ | sites of $C_5$, and $c$ | charges of $C_5$, and $2^{-44}$ |
| $C_7$ | $T$, $c\pm\tfrac14(1,1,1)$, $c\pm\tfrac1{2048}(1,1,1)$ | $1,1,1$, $q_5,q_5$, $q_7,q_7$ with $q_7=134713/2^{60}$ |

$C_5$ is the ABK construction: in units of the triangle's circumradius $\sqrt6$ the axial charges sit at
$\varepsilon = 1/(4\sqrt2)\approx 0.177$, and $q_5$ lies inside the window $[4316,4344]/2^{20}$ in which the
count is 24 (the ABK formula $\tfrac34\varepsilon^3-\tfrac5{32}\varepsilon^5$ gives $\approx 0.004118$).

$C_7$ repeats the construction at the centre of $C_5$ (ABK Remark 3 and Proposition 1). In circumradius
units the second pair sits at $\eta = 1/(2048\sqrt2)\approx 0.000345$. The centre of $C_5$ has
$\alpha \approx 0.00116$ and $\beta = 15/(8\sqrt6)$ in the expansion $U = U(0) + \alpha H_2 + \beta H_3 + O(4)$,
so ABK's choice $q_\eta = \alpha\eta^3 - \gamma\eta^5$ with $\gamma = \beta^2/(30\alpha)$ predicts
$q_7 \approx 1.168\times 10^{-13}$. A numerical scan at this $\eta$ finds the count 44 exactly for
$q_7 = k/2^{60}$ with $134615 \le k \le 134812$; $C_7$ uses the centre of that window, $k = 134713$.
At larger second-pair scales the outer trio of the bifurcated cluster does not appear (the counts found
at offsets $1/512$ and $1/1024$ per coordinate reach only 38), which matches the size of the degree-6
correction at the outer cluster radius.

### Theorem 1 (computer-certified)

| configuration | equilibria | index 1 | index 2 | all non-degenerate | $(n-1)^2$ |
|---|---|---|---|---|---|
| $C_3$ | 4 | 1 | 3 | yes | 4 |
| $C_4$ | **9** | 3 | 6 | yes | 9 |
| $C_5$ | **24** | 10 | 14 | yes | 16 |
| $C_6$ | **29** | 12 | 17 | yes | 25 |
| $C_7$ | **44** | 19 | 25 | yes | 36 |

"Equilibria" is the total number of critical points of $V$ in $\mathbb{R}^3$ minus the sites. Each certificate
encloses every critical point in a box that contains exactly one, proves that the Hessian is non-singular
there, and excludes the rest of space. The Morse relation $m_2-m_1=n-1$ holds in each case; it is not used
by the certificates, so it serves as a consistency check.

**Corollary.** Maxwell's bound is attained for four charges and exceeded for six and seven, in each case by
an explicit configuration with a completely counted, non-degenerate equilibrium set. It is now known to fail
for $n=5,6,7$.

### Proposition 2 (inserting a small charge; proof sketch)

Let $p$ be a non-degenerate critical point of index 1 of a potential $V$ of positive charges, with Hessian
eigenvalues $k_1\ge k_2>0>k_3$ and eigenvectors $e_1,e_2,e_3$. Add a charge $\kappa\lambda^3$ at $p$.

*(a) If $k_1>k_2$*, then for small $\lambda$ the point $p$ is replaced by exactly four non-degenerate critical
points near $\pm(\kappa/k_1)^{1/3}\lambda e_1$ (index 2) and $\pm(\kappa/k_2)^{1/3}\lambda e_2$ (index 1),
and all other critical points persist. The count rises by 3.

*(b) If $k_1=k_2$*, and the cubic term of $V$ at $p$, restricted to the circle in the $e_1e_2$-plane, is a
non-zero combination of $\cos3\theta$ and $\sin3\theta$ only (as under 3-fold rotational symmetry about
$e_3$), then $p$ is replaced by six non-degenerate critical points (three of each index). The count rises by 5.

*Sketch.* With $x=p+\lambda X$,
$V_{\text{new}}(x)=V(p)+\lambda^2 F(X)+\lambda^3R(\lambda,X)$ where $F(X)=\tfrac12X^{T}KX+\kappa/|X|$ and $R$
is smooth on compact subsets of $\mathbb{R}^3\setminus\{0\}$. $\nabla F=0$ forces $X$ to be an eigenvector of
$K$ with eigenvalue $\kappa/|X|^3>0$, which gives the four points in (a). At $X=s e_j$ the Hessian of $F$ is
$K-k_j(I-3e_je_j^{T})$, with eigenvalues $3k_j$ along $e_j$ and $k_i-k_j$ along $e_i$ ($i\ne j$). That gives
index 2 for $j=1$ and index 1 for $j=2$, both non-degenerate when $k_1\ne k_2$. The implicit function theorem
continues these to small $\lambda$. Close to $p$ the new charge's field dominates; between scale $\lambda$
and a fixed radius, $|\nabla V|\ge c|x-p|$ dominates. So no other critical points appear near $p$, and
critical points away from $p$ persist because the added potential tends to 0 in $C^2$ there. In case (b)
$F$ has a circle of critical points at $|X|=(\kappa/k)^{1/3}$ in the $e_1e_2$-plane. The first correction
along the circle is $\lambda^3$ times the cubic term, a trigonometric polynomial with six non-degenerate
critical points. The standard Morse–Bott perturbation argument then gives six critical points.

$C_4$ is case (b) applied to the centre of $C_3$ ($4-1+6=9$). $C_6$ is case (b) applied to the centre of
$C_5$, which by ABK's Lemma 2 is a $D_{3h}$-symmetric index-1 point ($24-1+6=29$). Theorem 1 does not depend on
Proposition 2; the certificates check the outcome directly.

**Consequences.** Combined with ABK's Proposition 1 (which keeps a $D_{3h}$ index-1 point at the origin after
every step), the maximum number $M(n)$ of non-degenerate equilibria satisfies

$$M(3+2m)\ge 4+20m,\qquad M(4+2m)\ge 9+20m .$$

These exceed $(n-1)^2$ only for $n=5,6,7$. They tie with it at $n=8$ ($49$) and $n=9$ ($64$), and fall short
for $n\ge10$. Every construction known to us is linear in $n$, while Maxwell's bound is quadratic. Whether
the conjecture fails for large $n$ is therefore open, as are $n=3$ and $n=4$.

## 3. Evidence for three and four charges

### Three charges: the inverse map

For a triangle $A_1A_2A_3$ and a point $x$ inside it there is exactly one charge vector, up to scale, for which
$x$ is an equilibrium. It follows from Lami's theorem:
$q_i=r_i^2 w_i$ with $(w_1:w_2:w_3)=(u_2\times u_3 : u_3\times u_1 : u_1\times u_2)$, where
$u_i=(x-A_i)/|x-A_i|$. Let $Z(x)$ be that charge ratio and $D(x)$ the determinant of the in-plane Hessian of
the corresponding potential at $x$. Differentiating $\nabla V(x;Z(x))=0$ shows that the folds of $Z$ are
exactly the degenerate equilibria, and that $Z$ preserves orientation on $R_{\min}=\{D>0\}$ (local minima)
and reverses it on saddles. So for each triangle:

> Maxwell's bound for three charges on this triangle $\iff$ $Z$ is one-to-one on $R_{\min}$.

`explore/three_scan.py` checks this over 1,289 triangle shapes in the fundamental domain (longest side 1,
heights from $10^{-4}$ to $0.87$). For each shape it rasterises $Z(R_{\min})$ onto a fine grid of the charge
simplex and counts how often each pixel is covered.

* In the 799 shapes where $R_{\min}$ is resolved, no charge ratio was covered twice, apart from three edge
  pixels. `explore/three_scan_check.py` re-solved those three directly and found 0 or 1 minimum, never 2;
  they are rasterisation artifacts at cusps of the fold image. Ninety sampled singly covered pixels all had
  exactly one minimum.
* In the other 487 shapes, mostly very flat triangles, no minimum region showed up at the grid's resolution.
  For those the scan is weak evidence: a grid cannot see a sliver thinner than its resolution.

Details are in `results/three_scan_summary.txt`.

Random sampling of 16,049 three-charge configurations (log-uniform charge ratios over 8 decades) found 2
with 4 equilibria and none with 6.

### Four charges

No configuration with more than 9 equilibria was found by any of these:

* the $C_{3v}$ family (triangle plus one axial charge), scanned finely around the value that cancels the
  quadratic term: maximum 9 (`explore/scan_c3v_fine.py`, 2,928 configurations);
* random multi-scale sampling of all 8 shape and charge parameters, about 25,000 configurations: maximum 7;
* ghost-guided hill climbing (8 chains, 300 steps), which tries to push near-miss minima of $|E|$ to zero:
  maximum 9;
* the four-charge inverse map, where each interior point of a tetrahedron fixes the charges that make it an
  equilibrium. Maxwell's bound 9 is equivalent to $m_1\le3$, so we searched by least squares for four
  distinct index-1 points sharing one charge vector (17 unknowns, 9 equations, 160 restarts). The residual
  never reached zero; the best runs flatten the tetrahedron and pull the points together, suggesting that
  the infimum is approached only in degenerate limits;
* importance sampling through the same inverse map (`explore/sample4_inverse.py`): 2,703 random tetrahedra
  and 16,179 charge vectors, each with at least one 1-saddle. The counts were 5, 7 or 9, never more.

Two local models explain why the natural designs stop where they do:

* **One small charge in a smooth background** (charge plus linear plus quadratic field) has at most 4
  equilibria near it in 3,000 random trials, or a ring broken into 6 when the background is axially
  symmetric (`explore/local_model.py`). So "three far charges plus one near charge" gives at most
  $3+6=9$ at leading order.
* **The ABK local model** $\Phi=L+Q+C+Z_4$ (harmonic linear, quadratic and cubic terms plus the axial pair's
  zonal quartic $Z_4=3\rho^4-24\rho^2z^2+8z^4$). The gradient is a cubic system, so Bézout allows 27
  solutions. But $\nabla Z_4$ vanishes on the complex isotropic directions $z=0$, $x^2+y^2=0$, so some
  solutions sit at infinity. Numerically a generic choice has 21 finite complex solutions (143 of 150
  random trials; the rest were solver misses), and at ABK's coefficients all 21 are real
  (`explore/local_poly.py`). In this sense ABK's construction is already optimal for "triangle plus axial
  pair", and 24 is the ceiling of that design.

## 4. How the certificates work and how to check them

Each program proves the same statement for a configuration: inside the padded bounding box of the sites
(which contains the convex hull, hence every equilibrium) $\nabla V$ has exactly $N$ zeros, each
non-degenerate, with the stated indices. For $C_3$ to $C_6$ both programs prove it. For $C_7$ only the
independent verifier currently finishes (see below), so that count rests on one program for now.

**Primary: `certify.py`**, using vectorised float64 interval arithmetic with every operation rounded
outward by one ulp (`ia.py`; IEEE operations are correctly rounded, so this is rigorous).
`tests/test_ia.py` checks `ia.py` against exact rational arithmetic. The box is bisected; a box is discarded
when

* the natural or the mean-value enclosure $g(c)+J(X)(X-c)$ of $\nabla V$ misses 0 in some component, or
* one site dominates: $q_i/r_{\max,i}^2>\sum_{j\ne i}q_j/r_{\min,j}^2$ on the box. This also covers boxes
  that contain a site.

The remaining boxes are inflated by 10% and passed to the Krawczyk operator
$K=c-Cg(c)+(I-CJ(X))(X-c)$. The test $K\subset\operatorname{int}X$ together with $\|I-CJ(X)\|_\infty<1$
proves a unique zero in $X$ at which every Hessian in $J(X)$ is non-singular. Overlapping verified boxes are
matched through tight Krawczyk enclosures. Morse indices come from Sylvester's law of inertia applied to
$Q^{T}J Q$ with Gershgorin discs.

**Independent: `verify_arb.py`.** It shares no code with the primary program: Arb ball arithmetic
(python-flint) for the field and Hessian, and exact rationals for box geometry and the dominance test. It
reads only the configuration and the list of zero boxes from the certificate, then:

1. re-proves each box with its own Krawczyk test and norm bound;
2. computes its own tight enclosures and checks that they are pairwise disjoint, so the zeros are distinct;
3. re-derives each Morse index;
4. runs its own bisection of the bounding box, discarding boxes by its own exclusion tests or by
   containment in a verified zero box. This proves that no other zeros exist.

```
pip install -r requirements.txt
python3 make_certs.py                      # writes the C3-C6 certificates  (about 5 s)
python3 explore/second_pair.py             # writes the C7 certificate      (seconds)
python3 verify_arb.py certificates/*.json  # independent check of all five  (about 2 min)
python3 tests/test_ia.py                   # interval arithmetic vs exact rationals
```

| certificate | primary: boxes, Krawczyk boxes | independent: boxes |
|---|---|---|
| $C_3$ | 2,759, 4 | 3,563 |
| $C_4$ | 8,311, 9 | 15,379 |
| $C_5$ | 92,429, 27 | 104,271 |
| $C_6$ | 112,203, 36 | 152,925 |
| $C_7$ | does not finish (issue #12) | 394,265 |

The $C_7$ boxes come from `explore/second_pair.py`: targeted Newton iteration seeded from the $C_5$
certificate points and the local-model cluster, then a Krawczyk box search and tightening in Arb
arithmetic. Nothing in that search needs to be trusted: `verify_arb.py` starts from the configuration
and the boxes alone and proves existence, uniqueness, indices and completeness as above. The float
certifier exceeds its box budget on $C_7$ (5,000,000 boxes in 37 s) and a raised-budget attempt was
killed at about 53 GB peak memory; making it finish is tracked in issue #12.

## 5. What is new and what is not

* The ABK construction, Morse counting, Krawczyk and Arb methods are prior work, credited above.
* $C_4$ is elementary. We could not find an earlier configuration attaining Maxwell's bound for $n=4$, and
  EFO state that none was known, but specialists may well have one.
* $C_5$ is ABK's example made explicit and counted completely.
* $C_6$ and Proposition 2 are simple combinations of known ideas. We found no earlier six-charge example
  above 25.
* $C_7$ is ABK's Proposition 1 at $m=2$ made explicit and counted completely. Its certificate is proved by
  the independent verifier only; the float certifier does not yet finish on it (issue #12).
* Novelty was checked by web search on 18 October 2026 only. Please tell us about earlier work.

**Who did what.** The question was chosen and all code, computation and writing for $C_3$ to $C_6$ were done
by Claude (an AI model by Anthropic) in one session on 18 October 2026, at the request of Foundation0
founding member Ahraz Arifuddin. The $C_7$ construction, search and certificate were added later the same
day by a second AI agent (amp-vela, Muse Code, muse-wonder-full-1.0). No person has reviewed the work
yet. The certificates are designed so that a reader does not have to trust the authors, only the verifiers,
which are short enough to read.

## 6. Next steps

* $n=4$: prove $M(4)=9$, or find a counterexample. The inverse map makes this a question about a map from a
  tetrahedron to the charge simplex.
* $n=3$: a computer-assisted proof on compact parts of the shape and charge space, through the
  inverse-map reformulation. The degenerate limits need separate asymptotic analysis.
* Superlinear constructions, which are needed to break Maxwell's bound for large $n$.
* A Lean statement of the certified counts (see google-deepmind/formal-conjectures issue 4682 on
  formalising the conjecture).

## References

* J. C. Maxwell, *A Treatise on Electricity and Magnetism*, vol. 1, §113, 1873.
* A. Gabrielov, D. Novikov, B. Shapiro, Mystery of point charges, Proc. LMS 95 (2007) 443–472.
* Y.-L. Tsai, Special cases of three point charges, Nonlinearity 24 (2011) 3299–3321,
  doi:10.1088/0951-7715/24/12/002.
* Y.-L. Tsai, Maxwell's conjecture on three point charges with equal magnitudes, Physica D 309 (2015) 86–98,
  doi:10.1016/j.physd.2015.07.007.
* V. Zolotov, Upper bounds for the number of isolated critical points via the Thom–Milnor theorem, Anal. Math.
  Phys. 13 (2023) 81.
* H. Edelsbrunner, C. Fillmore, G. Oliveira, Counting equilibria of the electrostatic potential, Proc. LMS 132
  (2026) e70163; arXiv:2501.05315.
* P. Arathoon, G. Ball, M. D. Kvalheim, The Maxwell conjecture is false, arXiv:2607.27197 (July 2026).
* A. Gabrielov, D. Novikov, T. Novikov, B. Shapiro, From 12 to 6: sharpening the three-charge bound in
  Maxwell's problem, arXiv:2607.28785 (July 2026).
* G. Oliveira, Critical points of point charge potentials along lines, arXiv:2609.00336 (2026).
* A. Enciso, D. Peralta-Salas, The finiteness conjecture for equilibria of electric fields generated by point
  charges of one sign, arXiv:2609.03965 (2026).
* X. Duan, Counterexamples to GNS Conjecture 1.8, arXiv:2609.09637 (2026).
* R. Krawczyk, Newton-Algorithmen zur Bestimmung von Nullstellen mit Fehlerschranken, Computing 4 (1969);
  F. Johansson, Arb: efficient arbitrary-precision midpoint-radius interval arithmetic, IEEE Trans. Comput.
  66 (2017).

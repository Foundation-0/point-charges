"""Pure-Z/Q Krawczyk & Bareiss Sylvester-inertia certifier for point-charge equilibria.

Zero external dependencies, zero IEEE-754 floats, zero bare asserts.
Uses outward-rounded exact dyadic interval arithmetic (int + math.isqrt at 2^-160)
combined with exact fraction-free Bareiss determinants and explicit rational
Lagrange congruence diagonalization (S^T H S = D over Q).

Independent clean-room certification of the three four-charge configurations
(n4-eleven-short, n4-eleven-flat, n4-eleven-jet) from Foundation-0/point-charges#22
(@zero-aska, commit da0aaa4f44780e27eef1ca44d80f4d743592095a), plus an analytic n=2
control and an exact Law-2b 11 -> 9 saddle-node fold wall across width 1/10^9.
"""
from __future__ import annotations

from fractions import Fraction
import math
from typing import Any

PREC_BITS = 160
SCALE = 1 << PREC_BITS
PRECOND_BITS = 80
PRECOND_SCALE = 1 << PRECOND_BITS
MAX_RATIONAL_CHARS = 256
MAX_RATIONAL_BITS = 1024
MAX_FORMAT_BITS = 4096

SOURCE_PROVENANCE: dict[str, Any] = {
    "repository": "Foundation-0/point-charges",
    "pr_number": 22,
    "commit": "da0aaa4f44780e27eef1ca44d80f4d743592095a",
    "tree": "92121f912928589a0e9930631c09d52a6d22d01a",
    "source_certificate_sha256": {
        "n4-eleven-short": "921842862facd0a837e3400cf8dd1fffe6e57efdafa9ff1c043fd013ba33efcc",
        "n4-eleven-flat": "a17b38e01f23a8535b076ce06d604cba1554f80713ab0bc5beec8e4c66cf6ebe",
        "n4-eleven-jet": "0d9e335eeeda62d4302da95793cb09c77c0bd7cffc9e4ac359953ca5aeeb32e8",
    },
}

CONFIG_SPECS: dict[str, dict[str, Any]] = {
    "n2-analytic-control": {
        "sites": (("-1", "0", "0"), ("1", "0", "0")),
        "charges": ("1", "1"),
        "centers": (("0", "0", "0"),),
        "half_width": "1/1048576",
        "expected_count": 1,
        "expected_m1": 0,
        "expected_m2": 1,
    },
    "n4-eleven-short": {
        "sites": (
            ("-1759/5000", "-27373/100000", "-4033/20000"),
            ("-2821/50000", "1394/3125", "-5359/25000"),
            ("-137/50000", "141/20000", "16917/20000"),
            ("20431/50000", "-4351/25000", "-19749/100000"),
        ),
        "charges": ("3339/12500", "2801/10000", "1", "13137/50000"),
        "centers": (
            ("-82495682691732444799/4722366482869645213696", "-61953998670928389063/4722366482869645213696", "73696570702606782327/4722366482869645213696"),
            ("-25510137469928055537/4722366482869645213696", "9991393537990858455/4722366482869645213696", "-3031789191450910465/590295810358705651712"),
            ("-22990416637916273653/4722366482869645213696", "-16998570266089695801/4722366482869645213696", "20150041307261517283/4722366482869645213696"),
            ("-7251674760299903451/2361183241434822606848", "108792624247493521853/4722366482869645213696", "75839974744016008273/4722366482869645213696"),
            ("-1043481278413748481/1180591620717411303424", "29057313108563981345/4722366482869645213696", "10177216404349378379/2361183241434822606848"),
            ("-1001929992758494087/4722366482869645213696", "1016283071960008889/4722366482869645213696", "912705144237886197/590295810358705651712"),
            ("335468776739472707/2361183241434822606848", "-1287942833425239621/2361183241434822606848", "-370091827138712777/147573952589676412928"),
            ("977021486638462183/2361183241434822606848", "-3874599611142829063/1180591620717411303424", "-8075216157375472801/2361183241434822606848"),
            ("18644150234245902171/4722366482869645213696", "424159802703740395/147573952589676412928", "-10728773370998477989/2361183241434822606848"),
            ("13776957286547515757/2361183241434822606848", "-11496942643033745211/4722366482869645213696", "10755176532822988727/2361183241434822606848"),
            ("45921329036510110395/2361183241434822606848", "-38021429291511047603/4722366482869645213696", "18049860593464595967/1180591620717411303424"),
        ),
        "half_width": "1/1099511627776",
        "expected_count": 11,
        "expected_m1": 4,
        "expected_m2": 7,
    },
    "n4-eleven-flat": {
        "sites": (
            ("-17589337229137/50000000000000", "-13687953100171/50000000000000", "-1260190803117/6250000000000"),
            ("-225642277629/4000000000000", "44609571816631/100000000000000", "-21432722866669/100000000000000"),
            ("-276303604643/100000000000000", "353357166213/50000000000000", "4229132458863/5000000000000"),
            ("40859773346807/100000000000000", "-17404480112807/100000000000000", "-19751913069353/100000000000000"),
        ),
        "charges": ("1669492820757/6250000000000", "280104254533/1000000000000", "1", "26274615972673/100000000000000"),
        "centers": (
            ("-88603987617940436797/4722366482869645213696", "-33327106150047015705/2361183241434822606848", "39449228903757021009/2361183241434822606848"),
            ("-15587245104037166189/4722366482869645213696", "58602524289814575015/2361183241434822606848", "81361633282847994095/4722366482869645213696"),
            ("-3547261850130335267/2361183241434822606848", "1338768103626239195/2361183241434822606848", "-8325675929642603105/4722366482869645213696"),
            ("-3478841878906068113/4722366482869645213696", "-348842292669077795/590295810358705651712", "414853205987894335/1180591620717411303424"),
            ("-843491020932790565/2361183241434822606848", "10473343959476655701/4722366482869645213696", "5725309688362715041/4722366482869645213696"),
            ("-826773154849306179/4722366482869645213696", "-520834377733731895/2361183241434822606848", "499625355458162783/4722366482869645213696"),
            ("-15872584315908221/4722366482869645213696", "100557996243434121/2361183241434822606848", "-5002554290081040389/4722366482869645213696"),
            ("567942498163171945/4722366482869645213696", "-3067146810310474483/2361183241434822606848", "-3497750679598141097/2361183241434822606848"),
            ("1377825077183910049/1180591620717411303424", "4164074445028445091/4722366482869645213696", "-8183181774733933047/4722366482869645213696"),
            ("5779991181715381511/4722366482869645213696", "-1421344202860684857/2361183241434822606848", "799240993097700789/1180591620717411303424"),
            ("2975469702779383943/147573952589676412928", "-19741256099258794089/2361183241434822606848", "9319548604572473303/590295810358705651712"),
        ),
        "half_width": "1/1099511627776",
        "expected_count": 11,
        "expected_m1": 4,
        "expected_m2": 7,
    },
    "n4-eleven-jet": {
        "sites": (
            ("127337715853/125000000000", "508164207341/500000000000", "15953647983/15625000000"),
            ("469666822527/500000000000", "-964491758367/1000000000000", "-38409881883/40000000000"),
            ("-484828868133/500000000000", "944250040577/1000000000000", "-483764649059/500000000000"),
            ("-38208806403/40000000000", "-957336625609/1000000000000", "467251181817/500000000000"),
        ),
        "charges": ("1", "170213348547/200000000000", "172818007371/200000000000", "209598953761/250000000000"),
        "centers": (
            ("-515570097032602075/1180591620717411303424", "-3187616775471627/2361183241434822606848", "-8377115537511175/4722366482869645213696"),
            ("-62513922419507883/4722366482869645213696", "-31252872838035315/2361183241434822606848", "-91716842128734165/4722366482869645213696"),
            ("-566838070580785/147573952589676412928", "-2248884493819311/590295810358705651712", "-155667317283129641/2361183241434822606848"),
            ("-6251117104017719/4722366482869645213696", "-1796218085601177211/4722366482869645213696", "-3585484126268795/2361183241434822606848"),
            ("5836388338447883/2361183241434822606848", "504823092424444265/4722366482869645213696", "1426603491301179/590295810358705651712"),
            ("14170459257154199/4722366482869645213696", "2997733259609069/1180591620717411303424", "2327114105014308873/2361183241434822606848"),
            ("78395147661108853/4722366482869645213696", "19621212305513615/1180591620717411303424", "59726018644642349/4722366482869645213696"),
            ("34452764798006889/1180591620717411303424", "56733786081633123/2361183241434822606848", "25233744222351504283/2361183241434822606848"),
            ("72940511514310377/2361183241434822606848", "10464979402426732767/590295810358705651712", "187016499308420999/4722366482869645213696"),
            ("65658247576098123/590295810358705651712", "11281846783843319/4722366482869645213696", "706727929311045/295147905179352825856"),
            ("17538607674861387073/1180591620717411303424", "122300930545338129/4722366482869645213696", "95366020187991617/2361183241434822606848"),
        ),
        "half_width": "1/1099511627776",
        "expected_count": 11,
        "expected_m1": 4,
        "expected_m2": 7,
    },
}

BIFURCATION_SPEC: dict[str, Any] = {
    "base_config": "n4-eleven-short",
    "perturbed_charge_index": 2,
    "tau_lo": "4239/125000000",
    "tau_hi": "33913/1000000000",
    "colliding_indices": (6, 7),
    "fold_normal": ("19579/250000", "-1", "-79319/250000"),
    "annihilation_box": (
        ("197/1000000", "63/312500"),
        ("-1561/1250000", "-6191/5000000"),
        ("-10611/5000000", "-21161/10000000"),
    ),
    "centers_lo": (
        ("-683074522969816585975/37778931862957161709568", "-1026148835461196808917/75557863725914323419136", "608850682642802563213/37778931862957161709568"),
        ("-336217871490301940329/75557863725914323419136", "65004052745718977265/37778931862957161709568", "-79977026224920613319/18889465931478580854784"),
        ("-250728952885102196819/75557863725914323419136", "-45497282557000060669/18889465931478580854784", "207610654253938109299/75557863725914323419136"),
        ("-59805507993149755817/18889465931478580854784", "224481439241808596001/9444732965739290427392", "312243872395026983659/18889465931478580854784"),
        ("-22889056332810624503/37778931862957161709568", "148774452330758876799/37778931862957161709568", "101111809773083424899/37778931862957161709568"),
        ("-18858800522423345479/75557863725914323419136", "27100814576663853121/75557863725914323419136", "76324229051336637371/75557863725914323419136"),
        ("15039954727533910915/75557863725914323419136", "-23427193450171141305/18889465931478580854784", "-160039891313200834041/75557863725914323419136"),
        ("15078834782728291269/75557863725914323419136", "-23551305232354830195/18889465931478580854784", "-160197402253409207177/75557863725914323419136"),
        ("227193062563613259897/75557863725914323419136", "162290252572744763077/75557863725914323419136", "-133744873735151360687/37778931862957161709568"),
        ("313860239499364894497/75557863725914323419136", "-32769509778027916743/18889465931478580854784", "58562184014017957563/18889465931478580854784"),
        ("1525630080907425298757/75557863725914323419136", "-631510091664921978747/75557863725914323419136", "598147326636073935605/37778931862957161709568"),
    ),
    "half_width_lo": "1/4398046511104",
    "centers_hi": (
        ("-42692197494722502877/2361183241434822606848", "-64134362314928387553/4722366482869645213696", "38053200759609074855/2361183241434822606848"),
        ("-2626682811306879927/590295810358705651712", "4062721131800798793/2361183241434822606848", "-19994110015848948621/4722366482869645213696"),
        ("-3917580641280477081/1180591620717411303424", "-5687068241064511601/2361183241434822606848", "12975430385910690283/4722366482869645213696"),
        ("-14951389411635936817/4722366482869645213696", "112240815159368048615/4722366482869645213696", "78061029430517122143/4722366482869645213696"),
        ("-2861086307171568133/4722366482869645213696", "18596433741706158071/4722366482869645213696", "12638703476047921819/4722366482869645213696"),
        ("-1178681865883659129/4722366482869645213696", "1693841637253919803/4722366482869645213696", "4770180102549591913/4722366482869645213696"),
        ("7099704073214312905/2361183241434822606848", "10143017398659155655/4722366482869645213696", "-4179485111293318091/1180591620717411303424"),
        ("19616017502500376465/4722366482869645213696", "-8192273822951164601/4722366482869645213696", "14640327784441348105/4722366482869645213696"),
        ("95351976205560938135/4722366482869645213696", "-39469420404476737751/4722366482869645213696", "74768485956751805651/4722366482869645213696"),
    ),
    "half_width_hi": "1/1099511627776",
    "max_exclusion_leaves": 128,
}


class PointChargeCertError(ValueError):
    """Explicit fail-closed certification error (never stripped by python3 -O)."""

    def __init__(self, code: str, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(f"{code}: {message}")


def _is_ascii_digits(text: str) -> bool:
    return len(text) > 0 and text.isascii() and text.isdigit()


def parse_exact_rational(raw: Any) -> Fraction:
    """Strictly parse an exact int, Fraction, or canonical rational string ('p' or 'p/q')."""
    if type(raw) is int:
        if raw.bit_length() > MAX_RATIONAL_BITS:
            raise PointChargeCertError("RATIONAL_TOO_LARGE", "integer exceeds bit limit")
        return Fraction(raw, 1)
    if type(raw) is Fraction:
        if raw.numerator.bit_length() > MAX_RATIONAL_BITS or raw.denominator.bit_length() > MAX_RATIONAL_BITS:
            raise PointChargeCertError("RATIONAL_TOO_LARGE", "fraction exceeds bit limit")
        return raw
    if type(raw) is not str:
        raise PointChargeCertError("INVALID_RATIONAL_TYPE", "rational input must be int, Fraction, or str")
    if len(raw) == 0 or len(raw) > MAX_RATIONAL_CHARS or not raw.isascii():
        raise PointChargeCertError("INVALID_RATIONAL_STR", "rational string empty, oversized, or non-ASCII")
    parts = raw.split("/")
    if len(parts) == 1:
        num_s, den_s = parts[0], "1"
    elif len(parts) == 2:
        num_s, den_s = parts[0], parts[1]
    else:
        raise PointChargeCertError("INVALID_RATIONAL_STR", "rational string contains multiple '/'")
    neg = num_s.startswith("-")
    digits_num = num_s[1:] if neg else num_s
    if not _is_ascii_digits(digits_num) or not _is_ascii_digits(den_s):
        raise PointChargeCertError("INVALID_RATIONAL_STR", "rational components must be ASCII digits")
    num_val = -int(digits_num) if neg else int(digits_num)
    den_val = int(den_s)
    if den_val <= 0:
        raise PointChargeCertError("ZERO_DENOMINATOR", "rational denominator must be positive")
    return Fraction(num_val, den_val)


def format_rational(val: Fraction) -> str:
    """Format an exact Fraction as 'p' or 'p/q' with bounded integer bit-length."""
    if type(val) is not Fraction:
        raise PointChargeCertError("INVALID_RATIONAL_TYPE", "expected exact Fraction")
    if val.numerator.bit_length() > MAX_FORMAT_BITS or val.denominator.bit_length() > MAX_FORMAT_BITS:
        raise PointChargeCertError("RATIONAL_TOO_LARGE", "fraction exceeds formatting bit limit")
    if val.denominator == 1:
        return str(val.numerator)
    return f"{val.numerator}/{val.denominator}"


# ------------------------------------------------------------------ outward-rounded dyadic intervals (2^-160)

Interval = tuple[int, int]


def frac_to_iv(val: Fraction) -> Interval:
    """Exact outward-rounded dyadic enclosure [lo/2^P, hi/2^P] of val in Q."""
    p = val.numerator
    q = val.denominator
    lo = (p << PREC_BITS) // q
    hi = -(((-p) << PREC_BITS) // q)
    return (lo, hi)


def iv_add(a: Interval, b: Interval) -> Interval:
    return (a[0] + b[0], a[1] + b[1])


def iv_sub(a: Interval, b: Interval) -> Interval:
    return (a[0] - b[1], a[1] - b[0])


def iv_mul(a: Interval, b: Interval) -> Interval:
    p1 = a[0] * b[0]
    p2 = a[0] * b[1]
    p3 = a[1] * b[0]
    p4 = a[1] * b[1]
    lo = min(p1, p2, p3, p4) >> PREC_BITS
    hi = -((-max(p1, p2, p3, p4)) >> PREC_BITS)
    return (lo, hi)


def iv_sq(a: Interval) -> Interval:
    if a[0] <= 0 and a[1] >= 0:
        mag = max(-a[0], a[1])
        return (0, -((-(mag * mag)) >> PREC_BITS))
    s0 = a[0] * a[0]
    s1 = a[1] * a[1]
    if s0 <= s1:
        return (s0 >> PREC_BITS, -((-s1) >> PREC_BITS))
    return (s1 >> PREC_BITS, -((-s0) >> PREC_BITS))


def iv_abs_max(a: Interval) -> int:
    return max(-a[0], a[1])


def iv_inv_r3_r5(u: Interval) -> tuple[Interval, Interval]:
    """Outward-rounded enclosures of u^(-3/2) and u^(-5/2) for strictly positive u."""
    if u[0] <= 0:
        raise PointChargeCertError("SITE_OVERLAP", "distance-squared interval is not strictly positive")
    s_lo = math.isqrt(u[0] << PREC_BITS)
    hi_arg = u[1] << PREC_BITS
    r = math.isqrt(hi_arg)
    s_hi = r if r * r == hi_arg else r + 1
    num3 = 1 << (3 * PREC_BITS)
    d3_hi = u[1] * s_hi
    d3_lo = u[0] * s_lo
    inv3 = (num3 // d3_hi, (num3 + d3_lo - 1) // d3_lo)
    num5 = 1 << (4 * PREC_BITS)
    d5_hi = u[1] * u[1] * s_hi
    d5_lo = u[0] * u[0] * s_lo
    inv5 = (num5 // d5_hi, (num5 + d5_lo - 1) // d5_lo)
    return inv3, inv5


def iv_det3(h: list[list[Interval]]) -> Interval:
    """Direct outward-rounded 3x3 determinant enclosure."""
    m00 = iv_sub(iv_mul(h[1][1], h[2][2]), iv_mul(h[1][2], h[2][1]))
    m01 = iv_sub(iv_mul(h[1][0], h[2][2]), iv_mul(h[1][2], h[2][0]))
    m02 = iv_sub(iv_mul(h[1][0], h[2][1]), iv_mul(h[1][1], h[2][0]))
    return iv_add(iv_sub(iv_mul(h[0][0], m00), iv_mul(h[0][1], m01)), iv_mul(h[0][2], m02))


def eval_grad_hess_iv(
    sites_iv: list[list[Interval]],
    charges_iv: list[Interval],
    box_iv: list[Interval],
) -> tuple[list[Interval], list[list[Interval]]]:
    """Evaluate rigorous enclosures of grad V(X) and Hess V(X) for V(x) = sum q_j / ||x - a_j||."""
    grad: list[Interval] = [(0, 0), (0, 0), (0, 0)]
    hess: list[list[Interval]] = [[(0, 0), (0, 0), (0, 0)] for _ in range(3)]
    for site, charge in zip(sites_iv, charges_iv):
        dx = [iv_sub(box_iv[k], site[k]) for k in range(3)]
        u = iv_add(iv_add(iv_sq(dx[0]), iv_sq(dx[1])), iv_sq(dx[2]))
        inv3, inv5 = iv_inv_r3_r5(u)
        q_inv3 = iv_mul(charge, inv3)
        q_inv5 = iv_mul(charge, inv5)
        three_q_inv5 = (3 * q_inv5[0], 3 * q_inv5[1])
        for r in range(3):
            grad[r] = iv_sub(grad[r], iv_mul(q_inv3, dx[r]))
            for c in range(r, 3):
                term = iv_mul(three_q_inv5, iv_mul(dx[r], dx[c]))
                if r == c:
                    term = iv_sub(term, q_inv3)
                hess[r][c] = iv_add(hess[r][c], term)
                if c != r:
                    hess[c][r] = hess[r][c]
    return grad, hess


# ------------------------------------------------------------------ exact rational Bareiss & Sylvester inertia

def _require_seq(val: Any, code: str, label: str) -> list[Any] | tuple[Any, ...]:
    if type(val) not in (list, tuple):
        raise PointChargeCertError(code, f"{label} must be exact list or tuple")
    return val


def _require_exact_fraction(val: Any, code: str, label: str) -> Fraction:
    if type(val) is not Fraction:
        raise PointChargeCertError(code, f"{label} must be exact Fraction")
    return val


def _require_square_frac_matrix(matrix: Any, expected_n: int | None = None) -> list[list[Fraction]]:
    rows = _require_seq(matrix, "INVALID_MATRIX", "matrix")
    n = len(rows)
    if n == 0:
        raise PointChargeCertError("INVALID_MATRIX", "matrix must be non-empty")
    if expected_n is not None and n != expected_n:
        raise PointChargeCertError("INVALID_MATRIX", f"matrix must be {expected_n}x{expected_n}")
    out: list[list[Fraction]] = []
    for row in rows:
        r_seq = _require_seq(row, "INVALID_MATRIX", "matrix row")
        if len(r_seq) != n:
            raise PointChargeCertError("INVALID_MATRIX", "matrix must be square")
        out.append([_require_exact_fraction(entry, "INVALID_MATRIX", "matrix entry") for entry in r_seq])
    return out


def bareiss_det_exact(matrix: list[list[Fraction]]) -> Fraction:
    """Exact fraction-free integer Bareiss determinant for any non-empty square Q-matrix."""
    clean = _require_square_frac_matrix(matrix)
    n = len(clean)
    int_rows: list[list[int]] = []
    denom_prod = 1
    for row in clean:
        row_lcm = 1
        for entry in row:
            row_lcm = math.lcm(row_lcm, entry.denominator)
        denom_prod *= row_lcm
        int_rows.append([entry.numerator * (row_lcm // entry.denominator) for entry in row])
    sign = 1
    prev_pivot = 1
    for k in range(n - 1):
        piv_row = -1
        for i in range(k, n):
            if int_rows[i][k] != 0:
                piv_row = i
                break
        if piv_row == -1:
            return Fraction(0, 1)
        if piv_row != k:
            int_rows[k], int_rows[piv_row] = int_rows[piv_row], int_rows[k]
            sign = -sign
        pivot = int_rows[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                int_rows[i][j] = (int_rows[i][j] * pivot - int_rows[i][k] * int_rows[k][j]) // prev_pivot
        prev_pivot = pivot
    return Fraction(sign * int_rows[n - 1][n - 1], denom_prod)


def mat3_inv_exact(mat: list[list[Fraction]]) -> tuple[list[list[Fraction]], Fraction]:
    """Exact 3x3 inverse over Q via Bareiss determinant and exact cofactors."""
    clean = _require_square_frac_matrix(mat, expected_n=3)
    det = bareiss_det_exact(clean)
    if det == 0:
        raise PointChargeCertError("SINGULAR_MATRIX", "3x3 matrix is singular")
    c00 = clean[1][1] * clean[2][2] - clean[1][2] * clean[2][1]
    c01 = -(clean[1][0] * clean[2][2] - clean[1][2] * clean[2][0])
    c02 = clean[1][0] * clean[2][1] - clean[1][1] * clean[2][0]
    c10 = -(clean[0][1] * clean[2][2] - clean[0][2] * clean[2][1])
    c11 = clean[0][0] * clean[2][2] - clean[0][2] * clean[2][0]
    c12 = -(clean[0][0] * clean[2][1] - clean[0][1] * clean[2][0])
    c20 = clean[0][1] * clean[1][2] - clean[0][2] * clean[1][1]
    c21 = -(clean[0][0] * clean[1][2] - clean[0][2] * clean[1][0])
    c22 = clean[0][0] * clean[1][1] - clean[0][1] * clean[1][0]
    inv = [
        [c00 / det, c10 / det, c20 / det],
        [c01 / det, c11 / det, c21 / det],
        [c02 / det, c12 / det, c22 / det],
    ]
    return inv, det


def verify_congruence_exact(
    mat: list[list[Fraction]],
    cong: list[list[Fraction]],
    diag: list[Fraction],
) -> None:
    """Self-verify det(S) != 0 and S^T A S == diag(d) exactly over Q."""
    clean_mat = _require_square_frac_matrix(mat, expected_n=3)
    clean_cong = _require_square_frac_matrix(cong, expected_n=3)
    d_seq = _require_seq(diag, "INVALID_MATRIX", "diag")
    if len(d_seq) != 3:
        raise PointChargeCertError("INVALID_MATRIX", "diag must have length 3")
    clean_diag = [_require_exact_fraction(v, "INVALID_MATRIX", "diag entry") for v in d_seq]
    if bareiss_det_exact(clean_cong) == 0:
        raise PointChargeCertError("CONGRUENCE_SINGULAR", "congruence transform is singular")
    for r in range(3):
        for c in range(3):
            acc = Fraction(0, 1)
            for p_i in range(3):
                for q_i in range(3):
                    acc += clean_cong[p_i][r] * clean_mat[p_i][q_i] * clean_cong[q_i][c]
            expected = clean_diag[r] if r == c else Fraction(0, 1)
            if acc != expected:
                raise PointChargeCertError("CONGRUENCE_MISMATCH", "S^T A S failed exact diagonal verification")


def sylvester_inertia_exact(mat: list[list[Fraction]]) -> tuple[int, int, int, list[Fraction]]:
    """Exact Lagrange symmetric congruence S^T A S = diag(d) over Q for a 3x3 symmetric matrix."""
    clean = _require_square_frac_matrix(mat, expected_n=3)
    for r in range(3):
        for c in range(r + 1, 3):
            if clean[r][c] != clean[c][r]:
                raise PointChargeCertError("ASYMMETRIC_MATRIX", "Sylvester inertia requires a symmetric matrix")
    work = [[clean[r][c] for c in range(3)] for r in range(3)]
    cong = [[Fraction(1 if r == c else 0, 1) for c in range(3)] for r in range(3)]
    diag: list[Fraction] = []
    for k in range(3):
        piv = -1
        for i in range(k, 3):
            if work[i][i] != 0:
                piv = i
                break
        if piv == -1:
            pair: tuple[int, int] | None = None
            for i in range(k, 3):
                for j in range(i + 1, 3):
                    if work[i][j] != 0:
                        pair = (i, j)
                        break
                if pair is not None:
                    break
            if pair is None:
                for _ in range(k, 3):
                    diag.append(Fraction(0, 1))
                break
            pi, pj = pair
            for idx in range(3):
                cong[idx][pi] += cong[idx][pj]
            for idx in range(k, 3):
                work[pi][idx] += work[pj][idx]
            for idx in range(k, 3):
                work[idx][pi] += work[idx][pj]
            piv = pi
        if piv != k:
            for idx in range(3):
                cong[idx][k], cong[idx][piv] = cong[idx][piv], cong[idx][k]
            work[k], work[piv] = work[piv], work[k]
            for idx in range(3):
                work[idx][k], work[idx][piv] = work[idx][piv], work[idx][k]
        d_k = work[k][k]
        diag.append(d_k)
        for i in range(k + 1, 3):
            factor = work[i][k] / d_k
            for idx in range(3):
                cong[idx][i] -= factor * cong[idx][k]
            for j in range(k + 1, 3):
                work[i][j] -= factor * work[k][j]
        for i in range(k + 1, 3):
            work[i][k] = Fraction(0, 1)
            work[k][i] = Fraction(0, 1)
    verify_congruence_exact(mat, cong, diag)
    n_pos = sum(1 for d in diag if d > 0)
    n_neg = sum(1 for d in diag if d < 0)
    n_zero = sum(1 for d in diag if d == 0)
    return n_pos, n_neg, n_zero, diag


# ------------------------------------------------------------------ Krawczyk & exclusion certification

def _validate_charges_and_sites(
    sites: list[list[Fraction]],
    charges: list[Fraction],
) -> tuple[list[list[Interval]], list[Interval]]:
    s_seq = _require_seq(sites, "INVALID_SYSTEM", "sites")
    c_seq = _require_seq(charges, "INVALID_SYSTEM", "charges")
    if len(s_seq) < 2 or len(c_seq) != len(s_seq):
        raise PointChargeCertError("INVALID_SYSTEM", "need at least 2 sites and matching charges")
    sites_iv: list[list[Interval]] = []
    for site in s_seq:
        pt = _require_seq(site, "INVALID_SYSTEM", "site")
        if len(pt) != 3:
            raise PointChargeCertError("INVALID_SYSTEM", "each charge site must be in R^3")
        sites_iv.append([frac_to_iv(_require_exact_fraction(coord, "INVALID_SYSTEM", "site coord")) for coord in pt])
    charges_iv: list[Interval] = []
    for q in c_seq:
        qf = _require_exact_fraction(q, "INVALID_SYSTEM", "charge")
        if qf <= 0:
            raise PointChargeCertError("NONPOSITIVE_CHARGE", "all point charges must be strictly positive")
        charges_iv.append(frac_to_iv(qf))
    return sites_iv, charges_iv


def _validate_box(box_frac: tuple[tuple[Fraction, Fraction], ...]) -> None:
    b_seq = _require_seq(box_frac, "INVALID_BOX", "box")
    if len(b_seq) != 3:
        raise PointChargeCertError("INVALID_BOX", "box must have 3 coordinate intervals")
    for pair in b_seq:
        p_seq = _require_seq(pair, "INVALID_BOX", "box interval")
        if len(p_seq) != 2:
            raise PointChargeCertError("INVALID_BOX", "box interval must have 2 endpoints")
        lo = _require_exact_fraction(p_seq[0], "INVALID_BOX", "box lo")
        hi = _require_exact_fraction(p_seq[1], "INVALID_BOX", "box hi")
        if lo >= hi:
            raise PointChargeCertError("INVALID_BOX", "box lower endpoint must be strictly less than upper endpoint")


def _verify_hessian_in_box(h_tilde_int: list[list[int]], h_box: list[list[Interval]]) -> None:
    for r in range(3):
        for c in range(3):
            if h_tilde_int[r][c] < h_box[r][c][0] or h_tilde_int[r][c] > h_box[r][c][1]:
                raise PointChargeCertError("HESSIAN_MIDPOINT_OUTSIDE", "symmetric midpoint Hessian outside H(X)")


def _check_krawczyk_contraction(rho_num: int) -> None:
    if rho_num >= SCALE:
        raise PointChargeCertError("KRAWCZYK_NOT_CONTRACTIVE", "contraction norm ||I - Y H(X)||_inf >= 1")


def _verify_harmonic_hessian(
    h_tilde: list[list[Fraction]],
    det_tilde: Fraction,
    h_box: list[list[Interval]],
) -> tuple[int, int, int, int, int]:
    n_pos, n_neg, n_zero, _ = sylvester_inertia_exact(h_tilde)
    if n_zero != 0 or (n_pos, n_neg) not in ((2, 1), (1, 2)):
        raise PointChargeCertError("EARNSHAW_VIOLATION", "non-degenerate harmonic equilibrium must have inertia (2,1,0) or (1,2,0)")
    det_box = iv_det3(h_box)
    if det_box[0] <= 0 and det_box[1] >= 0:
        raise PointChargeCertError("DET_INTERVAL_CONTAINS_ZERO", "interval Hessian determinant encloses zero")
    det_box_sign = 1 if det_box[0] > 0 else -1
    expected_sign = 1 if n_neg == 2 else -1
    det_tilde_sign = 1 if det_tilde > 0 else -1
    if det_box_sign != expected_sign or det_tilde_sign != expected_sign:
        raise PointChargeCertError("INERTIA_DET_MISMATCH", "Hessian determinant sign disagrees with Sylvester inertia")
    tr_box = iv_add(iv_add(h_box[0][0], h_box[1][1]), h_box[2][2])
    if tr_box[0] > 0 or tr_box[1] < 0:
        raise PointChargeCertError("HARMONIC_TRACE_MISMATCH", "interval trace does not enclose 0")
    return n_pos, n_neg, n_zero, det_box_sign, det_tilde_sign


def certify_equilibrium_box(
    sites: list[list[Fraction]],
    charges: list[Fraction],
    box_frac: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> dict[str, Any]:
    """Certify a unique, non-degenerate equilibrium inside box_frac via pure-Z/Q Krawczyk + Bareiss Sylvester."""
    sites_iv, charges_iv = _validate_charges_and_sites(sites, charges)
    _validate_box(box_frac)
    mid_frac = [(lo + hi) / 2 for lo, hi in box_frac]
    mid_iv = [frac_to_iv(v) for v in mid_frac]
    box_iv = [(frac_to_iv(lo)[0], frac_to_iv(hi)[1]) for lo, hi in box_frac]
    g_mid, h_mid = eval_grad_hess_iv(sites_iv, charges_iv, mid_iv)
    _, h_box = eval_grad_hess_iv(sites_iv, charges_iv, box_iv)

    h_tilde_int = [
        [((h_mid[r][c][0] + h_mid[r][c][1]) // 2 + (h_mid[c][r][0] + h_mid[c][r][1]) // 2) // 2 for c in range(3)]
        for r in range(3)
    ]
    _verify_hessian_in_box(h_tilde_int, h_box)
    h_tilde = [[Fraction(h_tilde_int[r][c], SCALE) for c in range(3)] for r in range(3)]
    y_exact, det_tilde = mat3_inv_exact(h_tilde)
    y_quant = [
        [Fraction((y_exact[r][c].numerator << PRECOND_BITS) // y_exact[r][c].denominator, PRECOND_SCALE) for c in range(3)]
        for r in range(3)
    ]
    y_iv = [[frac_to_iv(y_quant[r][c]) for c in range(3)] for r in range(3)]

    c_mat: list[list[Interval]] = [[(0, 0), (0, 0), (0, 0)] for _ in range(3)]
    rho_num = 0
    for r in range(3):
        row_sum = 0
        for c in range(3):
            acc = (0, 0)
            for k in range(3):
                acc = iv_add(acc, iv_mul(y_iv[r][k], h_box[k][c]))
            entry = iv_sub((SCALE, SCALE) if r == c else (0, 0), acc)
            c_mat[r][c] = entry
            row_sum += iv_abs_max(entry)
        if row_sum > rho_num:
            rho_num = row_sum
    _check_krawczyk_contraction(rho_num)

    dx_iv = [iv_sub(box_iv[k], mid_iv[k]) for k in range(3)]
    k_iv: list[Interval] = []
    for r in range(3):
        yg = (0, 0)
        cdx = (0, 0)
        for c in range(3):
            yg = iv_add(yg, iv_mul(y_iv[r][c], g_mid[c]))
            cdx = iv_add(cdx, iv_mul(c_mat[r][c], dx_iv[c]))
        k_iv.append(iv_add(iv_sub(mid_iv[r], yg), cdx))

    for d in range(3):
        k_lo = Fraction(k_iv[d][0], SCALE)
        k_hi = Fraction(k_iv[d][1], SCALE)
        if k_lo <= box_frac[d][0] or k_hi >= box_frac[d][1]:
            raise PointChargeCertError("KRAWCZYK_NOT_INTERIOR", "K(X) is not strictly contained in int(X)")

    n_pos, n_neg, n_zero, det_box_sign, det_tilde_sign = _verify_harmonic_hessian(h_tilde, det_tilde, h_box)

    return {
        "morse_index": n_neg,
        "inertia": [n_pos, n_neg, n_zero],
        "det_sign": det_box_sign,
        "krawczyk_rho": format_rational(Fraction(rho_num, SCALE)),
        "det_tilde_sign": det_tilde_sign,
        "harmonic_trace_contains_zero": True,
    }


def verify_pairwise_disjoint(
    boxes: list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]]],
) -> None:
    """Verify every pair of closed boxes is strictly disjoint in at least one coordinate."""
    b_seq = _require_seq(boxes, "INVALID_BOX", "boxes")
    for box in b_seq:
        _validate_box(box)
    count = len(b_seq)
    for i in range(count):
        for j in range(i + 1, count):
            separated = False
            for d in range(3):
                if b_seq[i][d][1] < b_seq[j][d][0] or b_seq[j][d][1] < b_seq[i][d][0]:
                    separated = True
                    break
            if not separated:
                raise PointChargeCertError("OVERLAPPING_BOXES", f"boxes {i} and {j} intersect")


def _box_strictly_inside(
    inner: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    outer: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> bool:
    for d in range(3):
        if inner[d][0] <= outer[d][0] or inner[d][1] >= outer[d][1]:
            return False
    return True


def _boxes_disjoint(
    a: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    b: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> bool:
    for d in range(3):
        if a[d][1] < b[d][0] or b[d][1] < a[d][0]:
            return True
    return False


def _centered_subbox_excluded(
    sites_iv: list[list[Interval]],
    charges_iv: list[Interval],
    subbox: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    v_iv: list[Interval],
) -> bool:
    mid_iv = [frac_to_iv((lo + hi) / 2) for lo, hi in subbox]
    box_iv = [(frac_to_iv(lo)[0], frac_to_iv(hi)[1]) for lo, hi in subbox]
    g_mid, _ = eval_grad_hess_iv(sites_iv, charges_iv, mid_iv)
    _, h_box = eval_grad_hess_iv(sites_iv, charges_iv, box_iv)
    dx_iv = [iv_sub(box_iv[k], mid_iv[k]) for k in range(3)]
    for r in range(3):
        comp = g_mid[r]
        for c in range(3):
            comp = iv_add(comp, iv_mul(h_box[r][c], dx_iv[c]))
        if comp[1] < 0 or comp[0] > 0:
            return True
    vg = (0, 0)
    for r in range(3):
        vg = iv_add(vg, iv_mul(v_iv[r], g_mid[r]))
    vh_dx = (0, 0)
    for c in range(3):
        vh_c = (0, 0)
        for r in range(3):
            vh_c = iv_add(vh_c, iv_mul(v_iv[r], h_box[r][c]))
        vh_dx = iv_add(vh_dx, iv_mul(vh_c, dx_iv[c]))
    proj = iv_add(vg, vh_dx)
    return proj[1] < 0 or proj[0] > 0


def certify_zero_free_box(
    sites: list[list[Fraction]],
    charges: list[Fraction],
    box_frac: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    fold_normal: list[Fraction],
    max_leaves: int = 128,
) -> dict[str, Any]:
    """Prove 0 not in grad V(x) for every x in box_frac via exact centered-form MVT bisection."""
    sites_iv, charges_iv = _validate_charges_and_sites(sites, charges)
    _validate_box(box_frac)
    fn_seq = _require_seq(fold_normal, "INVALID_NORMAL", "fold_normal")
    if len(fn_seq) != 3:
        raise PointChargeCertError("INVALID_NORMAL", "fold_normal must be a non-zero R^3 vector")
    fn_frac = [_require_exact_fraction(v, "INVALID_NORMAL", "fold_normal entry") for v in fn_seq]
    if all(v == 0 for v in fn_frac):
        raise PointChargeCertError("INVALID_NORMAL", "fold_normal must be a non-zero R^3 vector")
    if type(max_leaves) is not int or max_leaves < 1:
        raise PointChargeCertError("INVALID_BUDGET", "max_leaves must be a positive int")
    v_iv = [frac_to_iv(v) for v in fn_frac]
    stack = [box_frac]
    leaves = 0
    for _ in range(2 * max_leaves + 1):
        if not stack:
            return {"excluded": True, "leaves": leaves, "unresolved": 0}
        if leaves + len(stack) > max_leaves:
            return {"excluded": False, "leaves": leaves, "unresolved": len(stack)}
        current = stack.pop()
        if _centered_subbox_excluded(sites_iv, charges_iv, current, v_iv):
            leaves += 1
            continue
        widths = [current[d][1] - current[d][0] for d in range(3)]
        axis = 0
        if widths[1] > widths[axis]:
            axis = 1
        if widths[2] > widths[axis]:
            axis = 2
        mid = (current[axis][0] + current[axis][1]) / 2
        left_list = list(current)
        right_list = list(current)
        left_list[axis] = (current[axis][0], mid)
        right_list[axis] = (mid, current[axis][1])
        stack.append((left_list[0], left_list[1], left_list[2]))
        stack.append((right_list[0], right_list[1], right_list[2]))
    return {"excluded": False, "leaves": leaves, "unresolved": len(stack)}


def certify_configuration(
    name: str,
    sites: list[list[Fraction]],
    charges: list[Fraction],
    centers: list[list[Fraction]],
    half_width: Fraction,
    expected_count: int,
) -> tuple[dict[str, Any], list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]]]]:
    """Certify a full point-charge configuration and its exact Morse/Poincare-Hopf sum."""
    if type(name) is not str:
        raise PointChargeCertError("INVALID_NAME", "configuration name must be exact str")
    hw = _require_exact_fraction(half_width, "INVALID_HALF_WIDTH", "half_width")
    if hw <= 0:
        raise PointChargeCertError("INVALID_HALF_WIDTH", "half_width must be strictly positive")
    if type(expected_count) is not int or expected_count < 1:
        raise PointChargeCertError("COUNT_MISMATCH", "expected_count must be positive exact int")
    c_seq = _require_seq(centers, "INVALID_CENTER", "centers")
    if len(c_seq) != expected_count:
        raise PointChargeCertError("COUNT_MISMATCH", "center count disagrees with expected_count")
    boxes: list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]]] = []
    for ctr in c_seq:
        pt = _require_seq(ctr, "INVALID_CENTER", "center")
        if len(pt) != 3:
            raise PointChargeCertError("INVALID_CENTER", "each center must be in R^3")
        c0 = _require_exact_fraction(pt[0], "INVALID_CENTER", "center coord")
        c1 = _require_exact_fraction(pt[1], "INVALID_CENTER", "center coord")
        c2 = _require_exact_fraction(pt[2], "INVALID_CENTER", "center coord")
        boxes.append((
            (c0 - hw, c0 + hw),
            (c1 - hw, c1 + hw),
            (c2 - hw, c2 + hw),
        ))
    verify_pairwise_disjoint(boxes)
    m1 = 0
    m2 = 0
    max_rho = Fraction(0, 1)
    indices: list[int] = []
    for box in boxes:
        res = certify_equilibrium_box(sites, charges, box)
        idx = res["morse_index"]
        indices.append(idx)
        if idx == 1:
            m1 += 1
        else:
            m2 += 1
        rho = parse_exact_rational(res["krawczyk_rho"])
        if rho > max_rho:
            max_rho = rho
    poincare_hopf_sum = m2 - m1
    n_charges = len(sites)
    if poincare_hopf_sum != n_charges - 1:
        raise PointChargeCertError("MORSE_SUM_MISMATCH", "alternating Morse sum m2 - m1 != n - 1")
    maxwell_bound = (n_charges - 1) * (n_charges - 1)
    summary = {
        "name": name,
        "n_charges": n_charges,
        "maxwell_bound": maxwell_bound,
        "certified_equilibria_lower_bound": len(boxes),
        "exceeds_maxwell_bound": len(boxes) > maxwell_bound,
        "pairwise_disjoint": True,
        "all_nondegenerate": True,
        "index_counts": {"m0": 0, "m1": m1, "m2": m2, "m3": 0},
        "morse_alternating_sum": poincare_hopf_sum,
        "expected_poincare_hopf": n_charges - 1,
        "box_half_width": format_rational(half_width),
        "max_krawczyk_rho": format_rational(max_rho),
        "equilibrium_indices": indices,
    }
    return summary, boxes


def _load_spec_system(spec: dict[str, Any]) -> tuple[list[list[Fraction]], list[Fraction], list[list[Fraction]], Fraction]:
    sites = [[parse_exact_rational(v) for v in row] for row in spec["sites"]]
    charges = [parse_exact_rational(v) for v in spec["charges"]]
    centers = [[parse_exact_rational(v) for v in row] for row in spec["centers"]]
    half_width = parse_exact_rational(spec["half_width"])
    return sites, charges, centers, half_width


def _verify_bifurcation_checks(
    boxes_lo: list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]]],
    boxes_hi: list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]]],
    pre_indices: list[int],
    i_col: int,
    j_col: int,
    b_annih: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    excl: dict[str, Any],
) -> tuple[int, int]:
    idx_i = pre_indices[i_col]
    idx_j = pre_indices[j_col]
    if sorted((idx_i, idx_j)) != [1, 2]:
        raise PointChargeCertError("FOLD_INDEX_MISMATCH", "colliding pair must have opposite Morse indices (1 and 2)")
    if not _box_strictly_inside(boxes_lo[i_col], b_annih) or not _box_strictly_inside(boxes_lo[j_col], b_annih):
        raise PointChargeCertError("ANNIHILATION_BOX_MISMATCH", "colliding pre-fold boxes must lie strictly inside B_annih")
    for box in boxes_hi:
        if not _boxes_disjoint(box, b_annih):
            raise PointChargeCertError("SURVIVING_BOX_OVERLAP", "surviving post-fold box intersects B_annih")
    if not excl["excluded"]:
        raise PointChargeCertError("EXCLUSION_FAILED", "post-fold annihilation box was not certified zero-free")
    return idx_i, idx_j


def _check_expected_morse_counts(summary: dict[str, Any], expected_m1: int, expected_m2: int) -> None:
    if summary["index_counts"]["m1"] != expected_m1 or summary["index_counts"]["m2"] != expected_m2:
        raise PointChargeCertError("EXPECTED_INERTIA_MISMATCH", "unexpected Morse counts")


def certify_bifurcation_wall() -> dict[str, Any]:
    """Certify the Law-2b 11 -> 9 saddle-node fold wall across width 1/10^9 on n4-eleven-short."""
    base = CONFIG_SPECS[BIFURCATION_SPEC["base_config"]]
    sites = [[parse_exact_rational(v) for v in row] for row in base["sites"]]
    charges_base = [parse_exact_rational(v) for v in base["charges"]]
    q_idx = BIFURCATION_SPEC["perturbed_charge_index"]
    tau_lo = parse_exact_rational(BIFURCATION_SPEC["tau_lo"])
    tau_hi = parse_exact_rational(BIFURCATION_SPEC["tau_hi"])
    charges_lo = list(charges_base)
    charges_lo[q_idx] = charges_base[q_idx] + tau_lo
    charges_hi = list(charges_base)
    charges_hi[q_idx] = charges_base[q_idx] + tau_hi

    centers_lo = [[parse_exact_rational(v) for v in row] for row in BIFURCATION_SPEC["centers_lo"]]
    hw_lo = parse_exact_rational(BIFURCATION_SPEC["half_width_lo"])
    pre_summary, boxes_lo = certify_configuration("n4-short-pre-fold", sites, charges_lo, centers_lo, hw_lo, 11)

    centers_hi = [[parse_exact_rational(v) for v in row] for row in BIFURCATION_SPEC["centers_hi"]]
    hw_hi = parse_exact_rational(BIFURCATION_SPEC["half_width_hi"])
    post_summary, boxes_hi = certify_configuration("n4-short-post-fold", sites, charges_hi, centers_hi, hw_hi, 9)

    b_annih_raw = BIFURCATION_SPEC["annihilation_box"]
    b_annih = (
        (parse_exact_rational(b_annih_raw[0][0]), parse_exact_rational(b_annih_raw[0][1])),
        (parse_exact_rational(b_annih_raw[1][0]), parse_exact_rational(b_annih_raw[1][1])),
        (parse_exact_rational(b_annih_raw[2][0]), parse_exact_rational(b_annih_raw[2][1])),
    )
    i_col, j_col = BIFURCATION_SPEC["colliding_indices"]
    fold_normal = [parse_exact_rational(v) for v in BIFURCATION_SPEC["fold_normal"]]
    excl = certify_zero_free_box(sites, charges_hi, b_annih, fold_normal, BIFURCATION_SPEC["max_exclusion_leaves"])
    idx_i, idx_j = _verify_bifurcation_checks(
        boxes_lo, boxes_hi, pre_summary["equilibrium_indices"], i_col, j_col, b_annih, excl
    )

    colliding_linf = max(abs(centers_lo[i_col][d] - centers_lo[j_col][d]) for d in range(3))
    return {
        "base_config": BIFURCATION_SPEC["base_config"],
        "perturbed_charge_index": q_idx,
        "tau_lo": format_rational(tau_lo),
        "tau_hi": format_rational(tau_hi),
        "bracket_width": format_rational(tau_hi - tau_lo),
        "colliding_indices": [i_col, j_col],
        "colliding_morse_indices": [idx_i, idx_j],
        "colliding_linf_separation_at_tau_lo": format_rational(colliding_linf),
        "pre_fold": pre_summary,
        "post_fold": post_summary,
        "post_fold_zero_free_cover": {
            "annihilation_box": [[format_rational(lo), format_rational(hi)] for lo, hi in b_annih],
            "fold_normal": [format_rational(v) for v in fold_normal],
            "excluded": True,
            "leaves": excl["leaves"],
            "colliding_pair_strictly_inside": True,
            "surviving_nine_strictly_outside": True,
        },
        "poincare_hopf_invariant_across_fold": pre_summary["morse_alternating_sum"] == post_summary["morse_alternating_sum"],
    }


def build_full_payload() -> dict[str, Any]:
    """Compute the complete zero-float pure-Z/Q certification payload across all configurations."""
    configurations: dict[str, Any] = {}
    for name in ("n2-analytic-control", "n4-eleven-short", "n4-eleven-flat", "n4-eleven-jet"):
        spec = CONFIG_SPECS[name]
        sites, charges, centers, hw = _load_spec_system(spec)
        summary, _ = certify_configuration(name, sites, charges, centers, hw, spec["expected_count"])
        _check_expected_morse_counts(summary, spec["expected_m1"], spec["expected_m2"])
        configurations[name] = summary
    bifurcation = certify_bifurcation_wall()
    return {
        "track": "TRACK-S20A-POINT-CHARGES-N2-EQUILIBRIA",
        "arithmetic": {
            "mode": "pure_Z_Q_outward_dyadic_intervals_and_exact_Bareiss_Sylvester",
            "dyadic_precision_bits": PREC_BITS,
            "preconditioner_bits": PRECOND_BITS,
            "zero_floats": True,
        },
        "provenance": SOURCE_PROVENANCE,
        "configurations": configurations,
        "bifurcation_wall_11_to_9": bifurcation,
        "summary": {
            "n4_configurations_certified": 3,
            "total_n4_equilibria_certified": 33,
            "all_three_exceed_maxwell_bound_9": True,
            "fold_wall_bracket_width": bifurcation["bracket_width"],
            "scope_note": "Pure-Z/Q Krawczyk proves at least 11 pairwise-disjoint non-degenerate equilibria per configuration (sufficient to refute Maxwell's bound <= 9); global completeness exclusion over R^3 minus sites is out of scope of local Krawczyk.",
        },
    }

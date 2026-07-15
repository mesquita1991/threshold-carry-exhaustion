#!/usr/bin/env python3
"""Independent audit of the bridge B <-> R <-> A <-> D.

The script reconstructs the matrices directly from the definitions in the
manuscript.  Integer/rational parts are exact.  The confluent-root-of-unity
comparison is evaluated at 80 decimal digits and is used only to audit the
printed multiplicative factors and signs.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sp

x = sp.symbols("x")
mp.mp.dps = 80


def digit_coeff(q: int, n: int, u: int) -> int:
    if u < 0 or u > n * (q - 1):
        return 0
    coeffs = [1]
    for _ in range(n):
        nxt = [0] * (len(coeffs) + q - 1)
        for i, a in enumerate(coeffs):
            for d in range(q):
                nxt[i + d] += a
        coeffs = nxt
    return coeffs[u]


def scaled_B(q: int, n: int, t: int) -> sp.Matrix:
    N = q - 1
    return sp.Matrix(
        N,
        N,
        lambda i, j: digit_coeff(q, n, t - 2 - i - j)
        - digit_coeff(q, n, t - 1 - i - j),
    )


def exponent_data(q: int, n: int, t: int):
    N = q - 1
    m = n * N + 1
    h = t - N
    E = list(range(0, m - h)) + list(range(m - h + N, m + N))
    classes = {b: [e for e in E if e % q == b] for b in range(q)}
    return N, m, h, E, classes


def vandermonde(nodes: list[int]) -> int:
    return math.prod(y - x0 for x0, y in itertools.combinations(nodes, 2))


def polynomial_coefficients(poly: sp.Expr, n: int) -> list[sp.Expr]:
    P = sp.Poly(sp.expand(poly), x)
    return [P.coeff_monomial(x**k) for k in range(n)]


def construct_R_A(q: int, n: int, t: int):
    N, _, _, _, classes = exponent_data(q, n, t)
    k = {b: len(classes[b]) for b in range(q)}
    if any(k[b] > n for b in range(N)):
        return None
    d = {b: n - k[b] for b in range(N)}

    residual_columns: list[tuple[int, int, sp.Expr]] = []
    for b in range(N):
        w_b = sp.prod(x - e for e in classes[b])
        for j in range(d[b]):
            residual_columns.append((b, j, sp.expand(w_b * x**j)))

    E_N = classes[N]
    R = sp.Matrix(
        [
            [1]
            + [sp.expand(p.subs(x, 0) - p.subs(x, node)) for _, _, p in residual_columns]
            for node in E_N
        ]
    )

    A_columns: list[sp.Expr] = [sp.Integer(1)]
    A_columns.extend(p for _, _, p in residual_columns)
    w_N = sp.prod(x - e for e in E_N)
    for j in range(n - k[N]):
        A_columns.append(sp.expand(w_N * x**j))
    A = sp.Matrix(
        n,
        n,
        lambda i, j: polynomial_coefficients(A_columns[j], n)[i],
    )
    return R, A, classes, k, d


def pair_from_threshold(q: int, n: int, t: int):
    if t == q - 1:
        return n, n
    if t == (n + 1) * q - n:
        return 0, 1
    for j in range(2, n + 1):
        for s in range(n):
            if t == j * q - n + s:
                return n - j + 1, s + 1
    return None


def D_value(n: int, p: int, a0: int, q: int) -> sp.Integer:
    columns = []
    for a in range(1, n + 1):
        if a == a0:
            continue
        omitted = p - 1 if a < a0 else p
        W = sp.prod(x + a - j * q for j in range(n) if j != omitted)
        P = sp.Poly(sp.expand(W), x)
        columns.append([P.coeff_monomial(x**k) for k in range(1, n)])
    return sp.Integer(
        sp.Matrix(n - 1, n - 1, lambda i, j: columns[j][i]).det(method="domain-ge")
    )


def permutation_sign(indices: list[int]) -> int:
    inversions = sum(
        indices[i] > indices[j]
        for i in range(len(indices))
        for j in range(i + 1, len(indices))
    )
    return -1 if inversions % 2 else 1


def epsilons(q: int, n: int, E: list[int], classes, k, d):
    N = q - 1
    increasing = sorted(E)
    residue_order = [e for b in range(q) for e in sorted(classes[b])]
    eps_row = permutation_sign([increasing.index(e) for e in residue_order])
    exponent = N * (n - 1) + (k[N] - 1) + sum(
        d[b] * k[c] for b in range(N) for c in range(b + 1, N)
    )
    eps_col = -1 if exponent % 2 else 1
    return eps_row, eps_col


def complex_determinants(q: int, n: int, t: int):
    N, _, _, E, classes = exponent_data(q, n, t)
    zeta = mp.e ** (2j * mp.pi / q)
    C = []
    for e in E:
        row = [1]
        for j in range(1, N + 1):
            for s in range(n):
                row.append((e**s) * zeta ** (j * e))
        C.append(row)
    det_C = mp.det(mp.matrix(C))
    F = [[zeta ** (j * c) for j in range(q)] for c in range(q)]
    V = [[zeta ** (j * c) for j in range(1, q)] for c in range(N)]
    return det_C, mp.det(mp.matrix(F)), mp.det(mp.matrix(V)), E, classes


def confluent_denominator(q: int, n: int):
    N = q - 1
    zeta = mp.e ** (2j * mp.pi / q)
    y = [1] + [zeta**j for j in range(1, N + 1)]
    mu = [1] + [n] * N
    value = mp.mpf(1)
    value *= math.prod(math.factorial(s) for s in range(n)) ** N
    for j in range(1, N + 1):
        value *= y[j] ** (n * (n - 1) // 2)
    for a in range(N + 1):
        for b in range(a + 1, N + 1):
            value *= (y[b] - y[a]) ** (mu[a] * mu[b])
    return value


def main(output_path: Path) -> int:
    results: dict[str, object] = {}

    br_failures = []
    br_cases = 0
    ra_failures = []
    ra_cases = 0
    for q in range(2, 11):
        for n in range(2, 7):
            N = q - 1
            m = n * N + 1
            for t in range(N, N + m + 1):
                bdet = sp.Integer(scaled_B(q, n, t).det(method="domain-ge"))
                obj = construct_R_A(q, n, t)
                br_cases += 1
                if obj is None:
                    if bdet != 0:
                        br_failures.append(["capacity", q, n, t, int(bdet)])
                    continue
                R, A, classes, k, _ = obj
                rdet = sp.Integer(R.det(method="domain-ge"))
                if (bdet == 0) != (rdet == 0):
                    br_failures.append(["zero-equivalence", q, n, t, int(bdet), int(rdet)])
                ra_cases += 1
                rhs = (-1) ** (k[N] - 1) * vandermonde(classes[N]) * sp.Integer(
                    A.det(method="domain-ge")
                )
                if rdet != rhs:
                    ra_failures.append([q, n, t, int(rdet), int(rhs)])

    results["B_R_exact_zero_equivalence"] = {
        "cases": br_cases,
        "failures": br_failures,
    }
    results["R_A_exact_determinant_identity"] = {
        "cases": ra_cases,
        "failures": ra_failures,
    }

    ad_failures = []
    ad_cases = 0
    for n in range(2, 8):
        for q in range(n + 1, 13):
            thresholds = [q - 1]
            thresholds += [j * q - n + s for j in range(2, n + 1) for s in range(n)]
            thresholds += [(n + 1) * q - n]
            for t in thresholds:
                obj = construct_R_A(q, n, t)
                if obj is None:
                    ad_failures.append(["unexpected-capacity", q, n, t])
                    continue
                _, A, _, _, _ = obj
                pair = pair_from_threshold(q, n, t)
                if pair is None:
                    ad_failures.append(["missing-pair", q, n, t])
                    continue
                p, a0 = pair
                lhs = sp.Integer(A.det(method="domain-ge"))
                rhs = (-1) ** math.comb(n - 1, 2) * D_value(n, p, a0, q)
                ad_cases += 1
                if lhs != rhs:
                    ad_failures.append([q, n, t, p, a0, int(lhs), int(rhs)])

    results["A_D_exact_determinant_identity"] = {
        "cases": ad_cases,
        "failures": ad_failures,
    }

    bc_max_error = mp.mpf(0)
    bc_cases = 0
    cr_max_error = mp.mpf(0)
    cr_cases = 0
    printed_ce_r_counterexample = None
    printed_complete_counterexample = None

    for q in range(2, 8):
        for n in range(1, 4):
            N = q - 1
            m = n * N + 1
            if m > 13:
                continue
            for t in range(N, N + m + 1):
                det_C, det_F, det_V, E, classes = complex_determinants(q, n, t)
                C0 = confluent_denominator(q, n)
                h = t - N
                bdet = mp.mpf(int(scaled_B(q, n, t).det(method="domain-ge")))
                predicted_B = (-1) ** (N * (N + 1) // 2 + N * h) * det_C / C0
                bc_max_error = max(bc_max_error, abs(bdet - predicted_B))
                bc_cases += 1

                obj = construct_R_A(q, n, t)
                if obj is None:
                    continue
                R, _, _, k, d = obj
                det_R = mp.mpf(int(R.det(method="domain-ge")))
                if det_R == 0:
                    continue
                prod_delta = math.prod(vandermonde(classes[b]) for b in range(N))
                eps_row, eps_col = epsilons(q, n, E, classes, k, d)
                corrected = (
                    eps_row
                    * eps_col
                    * det_F
                    * det_V ** (n - 1)
                    * prod_delta
                    * det_R
                )
                cr_max_error = max(cr_max_error, abs(det_C / corrected - 1))
                cr_cases += 1

                printed_ce_r = (
                    eps_row
                    * eps_col
                    * det_F ** (-1)
                    * det_V ** (-(n - 1))
                    * prod_delta
                    * det_R
                )
                printed_complete = (
                    (-1) ** (N * (N + 1) // 2 + N * h)
                    * (q ** (-n * N))
                    * eps_row
                    * eps_col
                    * det_F
                    * det_V ** (n - 1)
                    / (C0 * prod_delta)
                    * det_R
                )
                true_unscaled_B = bdet * q ** (-n * N)
                if printed_ce_r_counterexample is None and abs(det_C - printed_ce_r) > mp.mpf("1e-30"):
                    printed_ce_r_counterexample = {
                        "q": q,
                        "n": n,
                        "t": t,
                        "printed_over_true": str(printed_ce_r / det_C),
                    }
                if printed_complete_counterexample is None and abs(true_unscaled_B - printed_complete) > mp.mpf("1e-30"):
                    printed_complete_counterexample = {
                        "q": q,
                        "n": n,
                        "t": t,
                        "true_det_B": str(true_unscaled_B),
                        "printed_rhs": str(printed_complete),
                        "printed_over_true": str(printed_complete / true_unscaled_B),
                    }

    results["B_C_confluent_identity_high_precision"] = {
        "cases": bc_cases,
        "maximum_absolute_error": str(bc_max_error),
    }
    results["C_R_corrected_factor_high_precision"] = {
        "cases": cr_cases,
        "maximum_ratio_error": str(cr_max_error),
        "corrected_formula": "det(C_E)=eps_row eps_col det(F+) det(V+)^(n-1) prod_{b<N}Delta_b det(R)",
    }
    results["printed_factor_counterexamples"] = {
        "CE_R_equation": printed_ce_r_counterexample,
        "complete_factor_equation": printed_complete_counterexample,
    }

    out = output_path
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(results, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=Path("audit_bridge_BRAD_results.json"),
        help="path for the deterministic JSON summary",
    )
    args = parser.parse_args()
    raise SystemExit(main(args.output))

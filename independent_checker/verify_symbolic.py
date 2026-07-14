#!/usr/bin/env python3
"""Exact symbolic regression checks for the main threshold-carry article.

This checker is deliberately independent of the C++ certificate generator.
It verifies only symbolic identities used in the paper; it does not attempt
to replace their proofs.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import sympy as sp

x, y, qsym = sp.symbols("x y q")


def D_polynomial(n: int, p: int, a0: int) -> sp.Poly:
    d = n - 1
    columns: list[list[sp.Expr]] = []
    for a in range(1, n + 1):
        if a == a0:
            continue
        omitted = p - 1 if a < a0 else p
        W = sp.prod(x + a - j * qsym for j in range(n) if j != omitted)
        poly = sp.Poly(sp.expand(W), x)
        columns.append([poly.coeff_monomial(x**k) for k in range(1, d + 1)])
    matrix = sp.Matrix(d, d, lambda i, j: columns[j][i])
    return sp.Poly(sp.expand(matrix.det(method="domain-ge")), qsym)


def D_numeric(n: int, p: int, a0: int, q_value: int) -> sp.Integer:
    d = n - 1
    columns: list[list[sp.Expr]] = []
    for a in range(1, n + 1):
        if a == a0:
            continue
        omitted = p - 1 if a < a0 else p
        W = sp.prod(x + a - j * q_value for j in range(n) if j != omitted)
        poly = sp.Poly(sp.expand(W), x)
        columns.append([poly.coeff_monomial(x**k) for k in range(1, d + 1)])
    return sp.Integer(sp.Matrix(d, d, lambda i, j: columns[j][i]).det(method="domain-ge"))


def expected_D0(n: int, a0: int) -> int:
    d = n - 1
    sign = -1 if (d * (d - 1) // 2) % 2 else 1
    binomial_product = math.prod(math.comb(d, k) for k in range(1, d + 1))
    vandermonde = math.prod(math.factorial(j) for j in range(1, n))
    return sign * binomial_product * vandermonde // (
        math.factorial(a0 - 1) * math.factorial(n - a0)
    )


def expected_D1_ratio(n: int, a0: int) -> sp.Rational:
    return sp.Rational(1, n - 1) * sum(
        sp.Rational(1, i - j)
        for i in range(a0 + 1, n + 1)
        for j in range(1, a0)
    )


def toeplitz_rhs_numeric(n: int, p: int, a0: int, q_value: int) -> sp.Integer:
    d = n - 1
    L = a0 - 1
    R = n - a0
    F = sp.expand(sp.prod(x - j * q_value for j in range(n) if j != p - 1))
    H = sp.expand(sp.prod(x - j * q_value for j in range(n) if j != p))

    basis: list[sp.Poly] = []
    current = F
    for _ in range(d + 1):
        basis.append(sp.Poly(current, x))
        current = sp.expand(current.subs(x, x + 1) - current)

    transform = sp.Matrix(
        [[basis[k].coeff_monomial(x**i) for k in range(d + 1)] for i in range(d + 1)]
    )
    target = sp.Matrix([sp.Poly(H, x).coeff_monomial(x**i) for i in range(d + 1)])
    connection = transform.inv() * target
    C = sum(connection[k] * y**k for k in range(d + 1))
    K = sp.Poly(sp.expand((1 + y) ** (L + 1) * C), y)

    def h(index: int) -> sp.Expr:
        if index < 0 or index > d:
            return 0
        return K.coeff_monomial(y**index)

    toeplitz = sp.Matrix(R, R, lambda i, j: h(L + i - j))
    gamma = math.prod(math.factorial(d) // math.factorial(d - k) for k in range(d))
    sign = -1 if (d * (d - 1) // 2) % 2 else 1
    return sp.Integer(sign * gamma) * toeplitz.det(method="domain-ge")


def main() -> int:
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    counts: dict[str, int] = {}
    cache: dict[tuple[int, int, int], sp.Poly] = {}

    for n in range(3, 7):
        for p in range(1, n):
            for a0 in range(2, n):
                poly = D_polynomial(n, p, a0)
                cache[(n, p, a0)] = poly
                assert int(poly.eval(0)) == expected_D0(n, a0)
                assert sp.Rational(poly.nth(1), poly.eval(0)) == expected_D1_ratio(n, a0)
                reflected = D_polynomial(n, n - p, n + 1 - a0)
                assert sp.expand(poly.as_expr() - reflected.as_expr()) == 0 or sp.expand(
                    poly.as_expr() + reflected.as_expr()
                ) == 0
                counts["stable_symbolic_families"] = counts.get("stable_symbolic_families", 0) + 1

    assert cache[(3, 1, 2)].as_expr() == -(qsym + 4)
    assert cache[(3, 2, 2)].as_expr() == -(qsym + 4)

    expected_n4 = {
        (1, 2): -5 * qsym**2 - 15 * qsym - 54,
        (1, 3): 3 * (qsym**2 - 5 * qsym - 18),
        (2, 2): -(qsym + 6) * (qsym + 9),
        (2, 3): -(qsym + 6) * (qsym + 9),
        (3, 2): 3 * (qsym**2 - 5 * qsym - 18),
        (3, 3): -5 * qsym**2 - 15 * qsym - 54,
    }
    for key, expected in expected_n4.items():
        assert sp.expand(cache[(4, *key)].as_expr() - expected) == 0

    exception = 48 * (qsym + 2) * (qsym**2 - 11 * qsym + 48)
    assert sp.expand(cache[(5, 1, 4)].as_expr() - exception) == 0
    assert sp.expand(cache[(5, 4, 2)].as_expr() - exception) == 0

    toeplitz_checks = 0
    for n in range(3, 8):
        for p in range(1, n):
            for a0 in range(2, n):
                for q_value in (2, 5):
                    assert D_numeric(n, p, a0, q_value) == toeplitz_rhs_numeric(
                        n, p, a0, q_value
                    )
                    toeplitz_checks += 1
    counts["toeplitz_exact_instances"] = toeplitz_checks
    counts["status"] = "passed"

    text = json.dumps(counts, indent=2, sort_keys=True) + "\n"
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

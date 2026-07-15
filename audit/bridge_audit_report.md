# Independent audit of the bridge `B <-> C_E <-> R <-> A <-> D`

## Status

Release `v1.0.1` corrects two multiplicative formulas printed in `v1.0.0` while preserving the vanishing equivalences used by the main theorem.

## Correct formulas

With the row and column conventions of the paper,

```text
det(C_E)
  = eps_row eps_col det(F+) det(V+)^(n-1)
    prod_{b<N} Delta_b det(R).
```

Consequently,

```text
det(B)
  = explicit nonzero Schur/confluent factor
    * eps_row eps_col det(F+) det(V+)^(n-1)
    * prod_{b<N} Delta_b / det(C_0)
    * det(R).
```

The `v1.0.0` intermediate formula left the Fourier determinants on the inverse side after solving for `det(C_E)`, and its complete formula placed the regular Vandermonde product in the denominator. Those equalities were false as printed.

## What remains unchanged

All factors outside `det(R)` are nonzero. Therefore

```text
det(B)=0 <=> det(C_E)=0 <=> det(R)=0 <=> det(A)=0,
```

and, in the stable regime,

```text
det(A)=0 <=> D_{n,p,a0}(q)=0.
```

The certified exhaustion theorem uses these vanishing equivalences and is unchanged.

## Audit artifacts

- `independent_checker/audit_bridge_BRAD.py`: direct reconstruction from definitions;
- `independent_checker/audit_bridge_BRAD_results.json`: frozen deterministic output;
- `audit/bridge_correction_patch.tex`: concise mathematical correction patch.

## Exact coverage

The script verifies:

- 990 `B`-to-`R` zero-equivalence cases;
- 615 exact `R`-to-`A` determinant identities;
- 790 exact stable `A`-to-`D` determinant identities.

It additionally performs 125 high-precision `B`-to-`C_E` checks and 64 high-precision corrected `C_E`-to-`R` factor checks. These complex numerical checks are regression tests for signs and factors, not substitutes for the symbolic proof in the paper.

The frozen output records no failures in any exact test and includes counterexamples to the two formulas printed in `v1.0.0`.

# Threshold Carry Events: corrected certified-exhaustion package

This repository contains the paper

> **Threshold Carry Events: Residual Determinants and Certified Exhaustion for Four to Thirty-Two Summands**

and the exact computational supplement used in its main theorem.

## Release v1.0.1 correction

Version `v1.0.1` corrects two multiplicative determinant formulas in the bridge from the confluent matrix `C_E` to the residual matrix `R`.

The correct identities are

```text
det(C_E)
  = eps_row eps_col det(F+) det(V+)^(n-1)
    prod_{b<N} Delta_b det(R),
```

and the corresponding complete `B`-to-`R` formula has `prod_{b<N} Delta_b` in the numerator. Version `v1.0.0` placed the Fourier determinants on the wrong side in the intermediate formula and placed the Vandermonde product in the denominator of the complete formula.

These printed equalities were false in `v1.0.0`. The singularity equivalence and the main exhaustion theorem are unchanged because all corrected factors are explicit and nonzero.

The paper now includes:

- a complete derivation of `B -> C_E` through a rectangular Schur specialization and the repeated-variable bialternant;
- an expanded `C_E -> R` Fourier/interpolation reduction with all determinant factors tracked;
- a table of spaces, bases, and dimensions;
- an explicit distinction between exact determinant identities and vanishing equivalences;
- an independent bridge-audit script, frozen JSON output, correction patch, and audit report.

## Main result

For independent uniform digits in base `q`, let `B_{q,r,t}` be the reduced matrix of the second-order Hoeffding interaction of the threshold event that the digit sum is at least `t`. The released proof establishes, for every base and threshold and for `4 <= r <= 32`, that singularity occurs exactly in the strict boundary ranges, the explicit carry bands, and the self-dual family.

The result is computer-assisted but exact. The symbolic part proves the chain

```text
Hoeffding matrix B
  -> confluent bialternant matrix C_E
  -> residual interpolation matrix R
  -> coefficient map A
  -> integer polynomial determinant D_{n,p,a0}(q).
```

The finite obligations are discharged by frozen finite-field certificates and an independent direct audit. Floating point is not used in the theorem-producing certificate checks. High-precision complex arithmetic appears only in the bridge regression audit for printed signs and normalization factors.

## One-command reproduction

Requirements:

- GNU/Linux or compatible shell;
- C++17 compiler with OpenMP;
- Python 3.11+;
- SymPy 1.14 and mpmath 1.3;
- `pdflatex` and an AMS LaTeX installation;
- Poppler utilities (`pdftotext`, `pdfinfo`);
- `sha256sum`.

Run:

```bash
bash ./run_all.sh
```

The command:

1. rebuilds the paper and verifies its 24-page structure and invariant text markers;
2. rebuilds all executables;
3. regenerates the stable certificate CSV and compares it byte-for-byte with the frozen archive;
4. verifies every no-root certificate independently;
5. reconstructs and checks every unstable parameter by two matrix-entry routes;
6. runs exact symbolic regression checks and archive-coverage tests;
7. regenerates the computational table used in the paper;
8. reruns the independent `B <-> R <-> A <-> D` audit and compares its JSON output with the frozen result;
9. verifies `MANIFEST.sha256`.

### Bridge-audit coverage

The independent audit records:

- `990` exact `B`-to-`R` zero-equivalence cases;
- `615` exact `R`-to-`A` determinant identities;
- `790` exact stable `A`-to-`D` determinant identities;
- `125` high-precision `B`-to-`C_E` regression cases;
- `64` high-precision corrected `C_E`-to-`R` factor cases;
- zero failures in all exact tests.

The numerical confluent tests are not used to prove the general identities. They are regression checks for the printed normalization, while the paper supplies the symbolic proofs.

### PDF reproducibility note

PDF object streams, producer metadata, line breaking, and embedded font subsets can differ across TeX Live releases even when the mathematical content is unchanged. The verification command therefore requires successful compilation, the expected page count, and invariant rendered-text markers. The tagged release separately publishes the exact PDF asset and its byte-level SHA-256 checksum.

## Repository structure

```text
paper/                 LaTeX source and generated computational table
src/                   certificate generator and table generator
certificates/          frozen stable no-root certificate archive
independent_checker/   certificate checkers and the independent bridge audit
audit/                 correction patch and human-readable audit report
tests/                 archive-coverage and regression tests
MANIFEST.sha256        integrity manifest
CITATION.cff           citation metadata
README.md              this file
run_all.sh             complete reproduction command
```

## Exact coverage of the exhaustion theorem

- stable mixed families, `5 <= n <= 30`: **8,112**;
- no-root prime certificates: **8,110**;
- symbolic exceptional families: **2**;
- unstable nonconstant triples: **112,375**;
- predicted unstable singular triples: **8,225**;
- certified unstable nonsingular triples: **104,150**.

## Scientific limits

The repository does **not** claim:

- exhaustion for arbitrary `r`;
- a uniform proof of nonvanishing for all residual degrees;
- a classification of the unresolved deep self-dual corank;
- a parameter-only global formula for inertia.

For each fixed residual degree, the paper proves finite decidability. The uniform obstacle is the nonvanishing of a family of rectangular Toeplitz-Pade minors in arbitrary degree.

## Citation

Use the metadata in [`CITATION.cff`](CITATION.cff). The tagged release identifies the exact corrected source, certificate set, audit files, and PDF used by the paper.

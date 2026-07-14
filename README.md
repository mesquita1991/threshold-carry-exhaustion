# Threshold Carry Events: certified exhaustion package

This repository contains the paper

> **Threshold Carry Events: Residual Determinants and Certified Exhaustion for Four to Thirty-Two Summands**

and the exact computational supplement used in its main theorem.

## Main result

For independent uniform digits in base `q`, let `B_{q,r,t}` be the reduced matrix of the second-order Hoeffding interaction of the threshold event that the digit sum is at least `t`. The released proof establishes, for every base and threshold and for `4 <= r <= 32`, that singularity occurs exactly in the strict boundary ranges, the explicit carry bands, and the self-dual family.

The result is computer-assisted but exact. The symbolic part proves the chain

```text
Hoeffding matrix B
  -> confluent-Fourier matrix
  -> residual interpolation matrix R
  -> coefficient map A
  -> integer polynomial determinant D_{n,p,a0}(q).
```

The finite obligations are discharged by frozen finite-field certificates and an independent direct audit. Floating point is not used.

## One-command reproduction

Requirements:

- GNU/Linux or compatible shell;
- C++17 compiler with OpenMP;
- Python 3.11+;
- SymPy 1.14;
- `pdflatex` and an AMS LaTeX installation;
- `sha256sum`.

Run:

```bash
./run_all.sh
```

The command:

1. rebuilds all executables;
2. regenerates the stable certificate CSV and compares it byte-for-byte with the frozen archive;
3. verifies every no-root certificate independently;
4. reconstructs and checks every unstable parameter by two matrix-entry routes;
5. runs exact symbolic regression checks and archive-coverage tests;
6. regenerates the computational table used in the paper;
7. rebuilds the PDF reproducibly and verifies its frozen SHA-256 digest;
8. verifies `MANIFEST.sha256`.

## Repository structure

```text
paper/                 LaTeX source, generated table, and frozen PDF digest
src/                   certificate generator and table generator
certificates/          frozen stable no-root certificate archive
independent_checker/   independent finite-field, unstable, and symbolic checkers
tests/                 archive-coverage and regression tests
MANIFEST.sha256        integrity manifest
CITATION.cff            citation metadata
README.md               this file
run_all.sh              complete reproduction command
```

## Exact coverage

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

For each fixed residual degree, the paper proves finite decidability. The uniform obstacle is the nonvanishing of a family of rectangular Toeplitz-Padé minors in arbitrary degree.

## Citation

Use the metadata in [`CITATION.cff`](CITATION.cff). The tagged release identifies the exact source and certificate set used by the paper; its release assets contain the verified PDF.

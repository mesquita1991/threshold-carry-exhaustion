#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$ROOT/build"
BIN="$BUILD/bin"
LOGS="$BUILD/logs"
CERTS="$BUILD/certificates"
PAPER_BUILD="$BUILD/paper"
PAPER_TMP="$(mktemp -d /tmp/threshold-carry-paper.XXXXXX)"
trap 'rm -rf "$PAPER_TMP"' EXIT

rm -rf "$BUILD"
mkdir -p "$BIN" "$LOGS" "$CERTS" "$PAPER_BUILD/generated" "$PAPER_TMP/generated"

for cmd in g++ python3 pdflatex pdftotext pdfinfo sha256sum cmp timeout; do
  command -v "$cmd" >/dev/null || { echo "missing required command: $cmd" >&2; exit 2; }
done


# Build the frozen paper before running the OpenMP verifiers.  Some TeX
# installations are more reliable when font loading occurs before the
# computational audit starts.
cp "$ROOT/paper/threshold_carry_exhaustion.tex" "$PAPER_TMP/"
cp "$ROOT/paper/generated/computational_counts.tex" "$PAPER_TMP/generated/"
(
  cd "$PAPER_TMP"
  export SOURCE_DATE_EPOCH=1784050033
  export FORCE_SOURCE_DATE=1
  : > "$PAPER_TMP/pdflatex.stdout.log"
  for pass in 1 2 3 4; do
    pdflatex -interaction=nonstopmode -halt-on-error threshold_carry_exhaustion.tex \
      >> "$PAPER_TMP/pdflatex.stdout.log" 2>&1
  done
)
cp "$PAPER_TMP/pdflatex.stdout.log" "$LOGS/pdflatex.log"
cp "$PAPER_TMP/threshold_carry_exhaustion.pdf" "$PAPER_BUILD/"

# PDF object streams and embedded font subsets vary across TeX Live releases.
# We therefore verify the invariant rendered content rather than requiring
# byte-identical PDF containers across installations.
pages="$(pdfinfo "$PAPER_BUILD/threshold_carry_exhaustion.pdf" \
  | awk '/^Pages:/ {print $2}')"
[[ "$pages" == "20" ]] || {
  echo "unexpected paper length: $pages pages (expected 20)" >&2
  exit 1
}
pdftotext -raw "$PAPER_BUILD/threshold_carry_exhaustion.pdf" \
  "$PAPER_BUILD/threshold_carry_exhaustion.txt"
for marker in \
  "THRESHOLD CARRY EVENTS:" \
  "Main Theorem." \
  "8,110" \
  "104,150" \
  "No claim is made for arbitrary r."; do
  grep -Fq "$marker" "$PAPER_BUILD/threshold_carry_exhaustion.txt" || {
    echo "missing expected rendered-text marker: $marker" >&2
    exit 1
  }
done

CXX="${CXX:-g++}"
CXXFLAGS=(-O3 -std=c++17)
OMPFLAGS=(-fopenmp)

"$CXX" "${CXXFLAGS[@]}" "${OMPFLAGS[@]}" \
  "$ROOT/src/generate_stable_certificates.cpp" \
  -o "$BIN/generate_stable_certificates"
"$CXX" "${CXXFLAGS[@]}" \
  "$ROOT/independent_checker/verify_stable_certificates.cpp" \
  -o "$BIN/verify_stable_certificates"
"$CXX" "${CXXFLAGS[@]}" "${OMPFLAGS[@]}" \
  "$ROOT/independent_checker/verify_unstable.cpp" \
  -o "$BIN/verify_unstable"

"$BIN/generate_stable_certificates" "$CERTS/regenerated.csv" \
  | tee "$LOGS/generate_stable_certificates.log"
cmp "$CERTS/regenerated.csv" \
  "$ROOT/certificates/stable_no_root_certificates_n5_n30.csv"

"$BIN/verify_stable_certificates" \
  "$ROOT/certificates/stable_no_root_certificates_n5_n30.csv" \
  | tee "$LOGS/verify_stable_certificates.log"

"$BIN/verify_unstable" "$LOGS/unstable_audit_summary.txt" \
  | tee "$LOGS/verify_unstable.log"

python3 "$ROOT/tests/test_small_cases.py" \
  | tee "$LOGS/tests.log"

python3 "$ROOT/src/generate_paper_table.py" \
  "$ROOT/certificates/stable_no_root_certificates_n5_n30.csv" \
  "$LOGS/unstable_audit_summary.txt" \
  "$PAPER_BUILD/generated/computational_counts.tex"
cmp "$PAPER_BUILD/generated/computational_counts.tex" \
  "$ROOT/paper/generated/computational_counts.tex"

# SymPy is deliberately last: this avoids leaving a large symbolic runtime
# ahead of the lightweight archive and table checks on constrained runners.
python3 "$ROOT/independent_checker/verify_symbolic.py" \
  "$LOGS/symbolic_audit_summary.json" \
  | tee "$LOGS/verify_symbolic.log"


(
  cd "$ROOT"
  sha256sum -c MANIFEST.sha256
)

echo "ALL CHECKS PASSED"

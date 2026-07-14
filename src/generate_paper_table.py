#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path


def parse_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def fmt(n: int) -> str:
    return f"{n:,}".replace(",", "{,}")


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: generate_paper_table.py CERT_CSV UNSTABLE_SUMMARY OUTPUT_TEX")
    cert_path, unstable_path, output_path = map(Path, sys.argv[1:])
    with cert_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    primes = [int(row["prime"]) for row in rows]
    unstable = parse_kv(unstable_path)

    per_n: dict[int, list[int]] = defaultdict(list)
    for row in rows:
        per_n[int(row["n"])].append(int(row["prime"]))

    summary = {
        "stable": len(rows),
        "certified": sum(p > 0 for p in primes),
        "exceptions": sum(p < 0 for p in primes),
        "maxprime": max(primes),
        "total": int(unstable["total"]),
        "predicted": int(unstable["predicted_singular"]),
        "nonpredicted": int(unstable["nonpredicted"]),
        "mismatch": int(unstable["mismatch"]),
        "route": int(unstable["route_mismatch"]),
    }

    lines = [
        r"\begin{center}",
        r"\small",
        r"\begin{tabular}{@{}lr@{}}",
        r"\toprule",
        r"Reproduced obligation & Count \\",
        r"\midrule",
        rf"Stable mixed families, $5\le n\le30$ & {fmt(summary['stable'])}\\",
        rf"No-root prime certificates & {fmt(summary['certified'])}\\",
        rf"Symbolic exceptional families & {fmt(summary['exceptions'])}\\",
        rf"Largest certificate prime & {fmt(summary['maxprime'])}\\",
        rf"Unstable nonconstant triples & {fmt(summary['total'])}\\",
        rf"Predicted unstable singular triples & {fmt(summary['predicted'])}\\",
        rf"Certified unstable nonsingular triples & {fmt(summary['nonpredicted'])}\\",
        rf"Classification mismatches & {fmt(summary['mismatch'])}\\",
        rf"Independent-entry mismatches & {fmt(summary['route'])}\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{center}",
        "",
        r"\begin{longtable}{@{}rrrr@{}}",
        r"\caption{Stable mixed-family certificate distribution by residual degree.}\label{tab:certificate-distribution}\\",
        r"\toprule",
        r"$n$ & Mixed families & Prime certificates & Largest prime\\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"$n$ & Mixed families & Prime certificates & Largest prime\\",
        r"\midrule",
        r"\endhead",
        r"\midrule",
        r"\multicolumn{4}{r}{Continued on next page}\\",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
    ]
    for n in sorted(per_n):
        vals = per_n[n]
        lines.append(
            rf"{n} & {fmt(len(vals))} & {fmt(sum(p > 0 for p in vals))} & {fmt(max(vals))}\\"
        )
    lines += [r"\end{longtable}", ""]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    expected = json.loads((ROOT / "tests" / "expected_counts.json").read_text())
    cert_path = ROOT / "certificates" / "stable_no_root_certificates_n5_n30.csv"
    with cert_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == expected["stable_mixed_families"]
    primes = [int(row["prime"]) for row in rows]
    assert sum(p > 0 for p in primes) == expected["stable_prime_certificates"]
    assert sum(p < 0 for p in primes) == expected["stable_symbolic_exceptions"]
    assert max(primes) == expected["maximum_certificate_prime"]

    exceptions = {
        (int(row["n"]), int(row["p"]), int(row["a0"]))
        for row in rows
        if int(row["prime"]) < 0
    }
    assert exceptions == {(5, 1, 4), (5, 4, 2)}

    seen = {(int(row["n"]), int(row["p"]), int(row["a0"])) for row in rows}
    expected_seen = {
        (n, p, a0)
        for n in range(5, 31)
        for p in range(1, n)
        for a0 in range(2, n)
    }
    assert seen == expected_seen
    print("small-case and archive-coverage tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Cross-package comparison benchmarks: OpenFermion (BosonOperator).

Scenario definitions, canonical keys, fairness contract and output format
live in BENCHMARKS.md.

OpenFermion notes (fairness contract, rule 5): coefficients must be numeric,
so this script uses ω = 1.0, J = 0.5, U = 0.25 where the other packages use
symbols. BosonOperator cannot be mixed with qubit/spin operators, so the
Jaynes–Cummings family, indexed sums, and symbolic expectation values are not
expressible; only the purely bosonic scenarios run.

Usage:
    python openfermion_bench.py      # needs `pip install openfermion`
"""

import datetime
import json
import os
import platform
import sys
import time
import timeit

import openfermion
from openfermion import BosonOperator, normal_ordered

CAP_SECONDS = 3.0
OMEGA, J, U = 1.0, 0.5, 0.25


def validate_equivalence():
    a = BosonOperator("0")
    ad = BosonOperator("0^")
    checks = [
        (
            "a a† = 1 + a†a",
            normal_ordered(a * ad) == BosonOperator("0^ 0") + BosonOperator(""),
        ),
        (
            "a a a† = a†a a + 2a",
            normal_ordered(a * a * ad)
            == BosonOperator("0^ 0 0") + 2 * BosonOperator("0"),
        ),
    ]
    ok = True
    for desc, passed in checks:
        ok &= passed
        print(f"  {'✓' if passed else '✗'}  {desc}", file=sys.stderr)
    return ok


def build_scenarios():
    scenarios = []

    # 2. Bosonic normal ordering (a·a†)ⁿ, eager left fold.
    a = BosonOperator("0")
    ad = BosonOperator("0^")

    def fock_reorder(n):
        acc = normal_ordered(a * ad)
        for _ in range(n - 1):
            acc = normal_ordered(acc * (a * ad))
        return acc

    for n in (2, 4, 6, 8, 10):
        scenarios.append((f"fock_reorder_n{n}", lambda n=n: fock_reorder(n)))

    # 3. Bose–Hubbard chain. Build-cost sweep over the chain length M feeds the
    # scaling figure (M = 8 reuses bose_hubbard_build); the fixed M = 8 pair
    # (build + H²) feeds the docs table.
    def bh_build_M(M):
        H = BosonOperator()
        for k in range(M):
            H += BosonOperator(f"{k}^ {k}", OMEGA)
        for k in range(M - 1):
            H += BosonOperator(f"{k}^ {k + 1}", J)
            H += BosonOperator(f"{k + 1}^ {k}", J)
        for k in range(M):
            H += BosonOperator(f"{k}^ {k}^ {k} {k}", U / 2)
        return normal_ordered(H)

    for M in (2, 4, 16):
        scenarios.append((f"bh_chain_M{M}", lambda M=M: bh_build_M(M)))

    Hbh = bh_build_M(8)
    scenarios.append(("bose_hubbard_build", lambda: bh_build_M(8)))
    scenarios.append(("bose_hubbard_H2", lambda: normal_ordered(Hbh * Hbh)))

    return scenarios


def time_scenario(thunk):
    thunk()  # warm-up
    t0 = time.perf_counter()
    thunk()
    trial = time.perf_counter() - t0
    if trial > CAP_SECONDS:
        return None, trial
    timer = timeit.Timer(thunk)
    loops, _ = timer.autorange()
    best = min(timer.repeat(repeat=5, number=loops)) / loops
    return best * 1e9, None


def main():
    print("Known-answer equivalence checks:", file=sys.stderr)
    if not validate_equivalence():
        raise SystemExit("equivalence checks failed; aborting")

    results, capped = {}, {}
    for key, thunk in build_scenarios():
        time_ns, cap = time_scenario(thunk)
        if time_ns is None:
            capped[key] = cap
            print(f"  {key}: capped ({cap:.1f} s/op)", file=sys.stderr)
        else:
            results[key] = {"time_ns": time_ns}
            print(f"  {key}: {time_ns / 1e6:.3f} ms", file=sys.stderr)

    out = {
        "meta": {
            "package": "OpenFermion",
            "version": openfermion.__version__,
            "language": "Python",
            "runtime_version": platform.python_version(),
            "date": datetime.date.today().isoformat(),
        },
        "results": results,
        "capped": capped,
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "openfermion.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {path}", file=sys.stderr)


if __name__ == "__main__":
    main()

"""Cross-package comparison benchmarks: OpenFermion (BosonOperator).

Scenario definitions, canonical keys, fairness contract and output format
live in BENCHMARKS.md.

OpenFermion notes (fairness contract, rule 5): coefficients must be numeric,
so this script uses ω = 1.0, J = 0.5, U = 0.25 where the other packages use
symbols. BosonOperator cannot be mixed with qubit/spin operators, so the
Jaynes–Cummings family, indexed sums, and symbolic expectation values are not
expressible; only the purely bosonic scenarios run.

Timing uses pyperf (a proper microbenchmark harness: per-process calibration,
warm-ups, and multiple worker processes); the reported figure is the minimum
per-call time over all samples, matching the other packages.

Usage:
    python openfermion_bench.py      # needs `pip install openfermion pyperf`
"""

import datetime
import json
import os
import platform
import sys
import tempfile
import time

import pyperf
import openfermion
from openfermion import BosonOperator, normal_ordered

CAP_SECONDS = 30.0
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

    for n in (2, 4, 6, 8, 10, 12, 14):
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

    for M in (2, 4, 16, 32):
        scenarios.append((f"bh_chain_M{M}", lambda M=M: bh_build_M(M)))

    Hbh = bh_build_M(8)
    scenarios.append(("bose_hubbard_build", lambda: bh_build_M(8)))
    scenarios.append(("bose_hubbard_H2", lambda: normal_ordered(Hbh * Hbh)))

    return scenarios


# Fair comparison: every harness gives each scenario the same ~10 s wall-clock
# budget and reports the minimum per-call time, in a single process. pyperf has
# no native time budget, so we size the number of values per scenario from a
# quick pre-timing so that (values x batch) ~ BUDGET_SECONDS, mirroring the Julia
# side's `run(b; seconds = 10)`. MAX_VALUES caps the sample count so fast
# scenarios stop early (~MAX_VALUES x MIN_TIME) instead of burning the full
# budget; only scenarios slower than BUDGET_SECONDS / MAX_VALUES per call use it
# all.
BUDGET_SECONDS = 10.0
MIN_TIME = 0.1  # pyperf per-value batch floor (seconds)
MAX_VALUES = 50  # early-stop sample cap for fast scenarios (~pyperf's default total)


def _values_for(op_time):
    return max(1, min(MAX_VALUES, round(BUDGET_SECONDS / max(op_time, MIN_TIME))))


def _single_trial(thunk):
    thunk()  # warm-up
    t0 = time.perf_counter()
    thunk()
    return time.perf_counter() - t0


def _capfile():
    return os.path.join(tempfile.gettempdir(), "sqa_capped_openfermion.json")


def main():
    runner = pyperf.Runner(processes=1, warmups=1, min_time=MIN_TIME)
    runner.parse_args()
    is_worker = bool(runner.args.worker)

    scenarios = build_scenarios()

    # The manager decides the >30 s cap and the per-scenario value count (both
    # from one quick timed call per scenario) and hands the capped set to
    # pyperf's worker subprocess through a temp file, so manager and worker
    # register the identical benchmark set. The value count reaches the worker
    # through pyperf's own `--values` argument, set per scenario below.
    if not is_worker:
        print("Known-answer equivalence checks:", file=sys.stderr)
        if not validate_equivalence():
            raise SystemExit("equivalence checks failed; aborting")
        capped, nvalues = {}, {}
        for key, thunk in scenarios:
            trial = _single_trial(thunk)
            if trial > CAP_SECONDS:
                capped[key] = trial
            else:
                nvalues[key] = _values_for(trial)
        with open(_capfile(), "w") as f:
            json.dump(list(capped), f)
    else:
        try:
            with open(_capfile()) as f:
                capped = {key: None for key in json.load(f)}
        except FileNotFoundError:
            capped = {}
        nvalues = {}

    results = {}
    for key, thunk in scenarios:
        if key in capped:
            continue
        if not is_worker:
            runner.args.values = nvalues[key]
        bench = runner.bench_func(key, thunk)
        if not is_worker and bench is not None:
            values = bench.get_values()
            if values:
                results[key] = {"time_ns": min(values) * 1e9}

    if is_worker:
        return

    try:
        os.remove(_capfile())
    except FileNotFoundError:
        pass

    for key in results:
        print(f"  {key}: {results[key]['time_ns'] / 1e6:.3f} ms", file=sys.stderr)
    for key, cap in capped.items():
        print(f"  {key}: capped ({cap:.1f} s/op)", file=sys.stderr)

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

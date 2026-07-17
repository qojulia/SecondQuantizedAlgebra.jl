"""Cross-package comparison benchmarks: SymPy (sympy.physics.quantum).

Scenario definitions, canonical keys, fairness contract and output format
live in BENCHMARKS.md.

SymPy note (fairness contract, rule 5): sympy's normal_ordered_form treats
boson and Pauli operators as mutually non-commuting, although they act on
different Hilbert spaces. `split_spaces` sorts boson factors before Pauli
factors in each product (a legal move since the spaces commute) so that like
terms merge; its cost is part of SymPy's timed canonicalization path.

Usage:
    python sympy_bench.py            # needs `pip install sympy`
"""

import datetime
import json
import platform
import sys
import time
import timeit

import sympy
from sympy import Mul, Pow, Rational, expand, symbols
from sympy.physics.quantum import Commutator, Dagger
from sympy.physics.quantum.boson import BosonOp
from sympy.physics.quantum.operatorordering import normal_ordered_form
from sympy.physics.quantum.pauli import (
    SigmaMinus,
    SigmaOpBase,
    SigmaPlus,
    qsimplify_pauli,
)

CAP_SECONDS = 3.0


def _is_pauli(f):
    base = f.base if isinstance(f, Pow) else f
    return isinstance(base, SigmaOpBase)


def split_spaces(expr):
    """Sort each product: scalars, then bosons, then Pauli operators."""
    expr = expr.expand()
    if expr.is_Add:
        return sum(split_spaces(t) for t in expr.args)
    c, nc = expr.args_cnc()
    bos = [f for f in nc if not _is_pauli(f)]
    pau = [f for f in nc if _is_pauli(f)]
    return Mul(*c) * Mul(*bos) * Mul(*pau)


def canon(expr):
    """Canonical form: commuting-spaces sort, boson normal order, Pauli algebra."""
    expr = split_spaces(expr)
    expr = normal_ordered_form(expr.expand(), independent=True, recursive_limit=20)
    return qsimplify_pauli(expr.expand()).expand()


def validate_equivalence():
    a = BosonOp("a")
    sp, sm = SigmaPlus(), SigmaMinus()
    checks = [
        ("a a† = 1 + a†a", expand(canon(a * Dagger(a)) - (Dagger(a) * a + 1)) == 0),
        (
            "a a a† = a†a a + 2a",
            expand(canon(a * a * Dagger(a)) - (Dagger(a) * a * a + 2 * a)) == 0,
        ),
        ("σ⁻σ⁺ + σ⁺σ⁻ = 1", expand(canon(sm * sp + sp * sm)) == 1),
    ]
    ok = True
    for desc, passed in checks:
        ok &= passed
        print(f"  {'✓' if passed else '✗'}  {desc}", file=sys.stderr)
    return ok


def build_scenarios():
    scenarios = []

    # 1. Jaynes–Cummings family.
    a = BosonOp("a")
    sp, sm = SigmaPlus(), SigmaMinus()
    wc, wa, g = symbols("omega_c omega_a g", real=True)

    def jc_build():
        return canon(wc * Dagger(a) * a + wa * sp * sm + g * (Dagger(a) * sm + a * sp))

    H = jc_build()
    scenarios.append(("jc_build", jc_build))
    scenarios.append(("jc_H2", lambda: canon(H * H)))
    scenarios.append(("jc_heisenberg", lambda: canon(Commutator(H, a).doit())))

    def nested(depth):
        x = sm
        for _ in range(depth):
            x = canon(Commutator(H, x).doit())
        return x

    for d in (1, 2, 3, 4, 5, 6):
        scenarios.append((f"jc_nested_d{d}", lambda d=d: nested(d)))

    def jc_pow(p):
        acc = H
        for _ in range(p - 1):
            acc = canon(acc * H)
        return acc

    for p in (2, 3, 4, 5, 6):
        scenarios.append((f"jc_pow_n{p}", lambda p=p: jc_pow(p)))

    # 2. Bosonic normal ordering (a·a†)ⁿ, eager left fold.
    def fock_reorder(n):
        acc = canon(a * Dagger(a))
        for _ in range(n - 1):
            acc = canon(acc * (a * Dagger(a)))
        return acc

    for n in (2, 4, 6, 8, 10):
        scenarios.append((f"fock_reorder_n{n}", lambda n=n: fock_reorder(n)))

    # 3. Bose–Hubbard chain. Build-cost sweep over the chain length M feeds the
    # scaling figure (M = 8 reuses bose_hubbard_build); the fixed M = 8 pair
    # (build + H²) feeds the docs table.
    w, J, U = symbols("omega J U", real=True)

    def bh_build_M(M):
        ops = [BosonOp(f"a{k}") for k in range(M)]
        H = sum(w * Dagger(ops[k]) * ops[k] for k in range(M))
        H += J * sum(
            Dagger(ops[k]) * ops[k + 1] + Dagger(ops[k + 1]) * ops[k]
            for k in range(M - 1)
        )
        H += (
            U
            * Rational(1, 2)
            * sum(Dagger(ops[k]) * Dagger(ops[k]) * ops[k] * ops[k] for k in range(M))
        )
        return canon(H)

    for M in (2, 4, 16):
        scenarios.append((f"bh_chain_M{M}", lambda M=M: bh_build_M(M)))

    Hbh = bh_build_M(8)
    scenarios.append(("bose_hubbard_build", lambda: bh_build_M(8)))
    scenarios.append(("bose_hubbard_H2", lambda: canon(Hbh * Hbh)))

    # 4./5. Indexed sums and symbolic expectation values: not expressible.
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
            "package": "SymPy",
            "version": sympy.__version__,
            "language": "Python",
            "runtime_version": platform.python_version(),
            "date": datetime.date.today().isoformat(),
        },
        "results": results,
        "capped": capped,
    }
    import os

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "sympy.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {path}", file=sys.stderr)


if __name__ == "__main__":
    main()

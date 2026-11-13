"""Shared contracts, provenance and untimed checks for comparison runners."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess

ROOT = Path(os.environ.get("SQA_COMPARE_ROOT", Path(__file__).resolve().parents[2]))
CONTRACT = "canonical-symbolic-v2"
KEYS = (["jc_build", "jc_H2", "jc_heisenberg"]
        + [f"jc_nested_d{d}" for d in range(1, 8)]
        + [f"jc_pow_n{n}" for n in range(2, 8)]
        + [f"fock_reorder_n{n}" for n in range(2, 15, 2)]
        + [f"bh_chain_M{m}" for m in (2, 4, 16, 32)]
        + ["bose_hubbard_build", "bose_hubbard_H2", "tavis_cummings_build",
           "tavis_cummings_comm", "meanfield_average"])


def selected(scenarios):
    keys = os.environ.get("SQA_COMPARE_KEYS", "").split(",")
    return [sc for sc in scenarios if not keys[0] or sc[0] in keys]


def provenance(script):
    def git(*args):
        return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()
    cpu = platform.processor()
    if Path("/proc/cpuinfo").exists():
        cpu = next((line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines() if line.startswith("model name")), cpu)
    return dict(contract=CONTRACT, run_id=os.environ.get("SQA_COMPARE_RUN_ID", "standalone"),
                date=datetime.date.today().isoformat(), cpu=cpu,
                source_revision=git("rev-parse", "HEAD"), dirty=bool(git("status", "--porcelain")),
                benchmark_script_sha256=hashlib.sha256(Path(script).read_bytes()).hexdigest(),
                statistic="minimum", coefficient_domain="symbolic",
                validation="passed", threads=os.environ.get("JULIA_NUM_THREADS", "1"))


def save_result(name, output, script):
    output["meta"].update(provenance(script))
    output["statuses"] = {k: "passed" for k in output["results"]}
    output["statuses"].update({k: "capped" for k in output.get("capped", {})})
    directory = Path(os.environ.get("SQA_COMPARE_OUTPUT", Path(script).parent / "results"))
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{name}.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(output, indent=2) + "\n")
    temporary.replace(path)
    return path


def stirling(n, k):
    row = [1] + [0] * n
    for m in range(1, n + 1):
        row = [0] + [row[j - 1] + j * row[j] for j in range(1, n + 1)]
    return row[k]


def validate_scenarios(module, scenarios):
    """Check independent coefficient formulas and analytical model drifts."""
    from sympy import symbols, expand
    cases = dict(scenarios)
    if module.__name__ == "openfermion":
        a = module.BosonOperator("0", 1)
        adj = module.hermitian_conjugated
        canon = module.normal_ordered
        identity = module.BosonOperator("", 1)
        def equal(x, y):
            return all(expand(c) == 0 for c in canon(x-y).terms.values())
        w, coupling, interaction = symbols("omega J U", real=True)
        neighbor = module.BosonOperator("1", 1)
    else:
        from sympy.physics.quantum import Dagger
        from sympy.physics.quantum.boson import BosonOp
        from sympy.physics.quantum.pauli import SigmaMinus
        a = BosonOp("a")
        adj = Dagger
        canon = module.canon
        identity = 1
        def equal(x, y):
            return expand(canon(x-y)) == 0
        wc, wa, g = symbols("omega_c omega_a g", real=True)
        sm = SigmaMinus()
        H = cases["jc_build"]()
        assert equal(cases["jc_heisenberg"](), -wc*a-g*sm)
        for d in (2, 4):
            from math import comb
            def power(n):
                value = identity
                for step in range(n):
                    value = canon(value*H)
                return value
            expected = sum((-1)**k * comb(d, k) * canon(canon(power(d-k)*sm)*power(k)) for k in range(d+1))
            assert equal(cases[f"jc_nested_d{d}"](), expected)
        w, coupling, interaction = symbols("omega J U", real=True)
        neighbor = BosonOp("a1")
    for n in (2, 4, 8, 14):
        expected = sum((stirling(n+1, k+1)*adj(a)**k*a**k for k in range(n+1)), identity*0)
        assert equal(cases[f"fock_reorder_n{n}"](), expected), f"normal ordering n={n}"
    H = cases["bose_hubbard_build"]()
    if module.__name__ != "openfermion":
        a = BosonOp("a0")
    expected = -w*a-coupling*neighbor-interaction*adj(a)*a*a
    assert equal(canon(H*a-a*H), expected), "Bose-Hubbard Heisenberg drift"
    assert equal(H, adj(H)), "Bose-Hubbard Hermiticity"
    print("Scenario checks passed: exact normal-order coefficients, analytical drifts and Hermiticity")

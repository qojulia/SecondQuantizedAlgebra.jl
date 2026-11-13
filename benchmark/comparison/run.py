"""Serial comparison runs with isolated output, provenance and publication checks."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import tarfile
from comparison_support import CONTRACT, KEYS, ROOT, provenance

HERE = Path(__file__).resolve().parent


def execute(command, environment, log, timeout):
    with log.open("w") as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            return "passed" if process.wait(timeout=timeout) == 0 else "failed"
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            return "timed_out"


def prepare(args):
    python = HERE / ".venv/bin/python"
    lock = HERE / "requirements.lock"
    if args.from_run:
        archived = args.from_run.resolve()
        snapshot = json.loads((archived / "run.json").read_text())
        revision = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
        if revision != snapshot["source_revision"]:
            raise ValueError("Check out the archived SQA revision before frozen reproduction")
        current_patch = subprocess.check_output(["git", "-C", str(ROOT), "diff", "HEAD", "--", "src", "ext", "Project.toml"], text=True)
        if current_patch != (archived / "source.patch").read_text():
            raise ValueError("Restore the archived SQA source changes before frozen reproduction")
        actual_julia = subprocess.check_output([args.julia, "--version"], text=True).strip()
        if snapshot.get("julia_runtime") and actual_julia != snapshot["julia_runtime"]:
            raise ValueError("Select the archived Julia runtime with --julia")
        if python.exists():
            actual_python = subprocess.check_output([str(python), "-c", "import platform; print(platform.python_version())"], text=True).strip()
            if actual_python != snapshot["python_runtime"]:
                raise ValueError("Recreate the comparison .venv with the archived Python runtime")
        for name in ("julia_bench.jl", "sympy_bench.py", "openfermion_bench.py", "sneg_bench.wls", "comparison_support.py"):
            if (HERE / name).read_bytes() != (archived / "source" / name).read_bytes():
                raise ValueError(f"Restore the archived comparison source before reproduction: {name}")
        lock = archived / "python-requirements.txt"
    if not python.exists():
        command = ["uv", "venv", str(python.parents[1])]
        if args.from_run:
            command += ["--python", snapshot["python_runtime"]]
        subprocess.run(command, check=True)
    if args.from_run:
        shutil.copy(archived / "Manifest.toml", ROOT / "benchmark/Manifest.toml")
    if args.refresh:
        subprocess.run(["uv", "pip", "compile", str(HERE / "requirements.in"), "--upgrade", "-o", str(lock)], check=True)
    subprocess.run(["uv", "pip", "sync", "--python", str(python), str(lock)], check=True)
    command = [args.julia, "--startup-file=no", "--project=benchmark", str(HERE / "prepare.jl")]
    if args.refresh:
        command.append("refresh")
    elif args.from_run:
        command.append("frozen")
    subprocess.run(command, cwd=ROOT, check=True)


def enrich(directory, packages, context):
    for name in packages:
        path = directory / f"{name}.json"
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        data["meta"].update(context)
        scripts = {"sqa": "julia_bench.jl", "quantumalgebra": "julia_bench.jl", "sympy": "sympy_bench.py", "openfermion": "openfermion_bench.py", "sneg": "sneg_bench.wls"}
        data["meta"]["benchmark_script_sha256"] = hashlib.sha256((directory / "source" / scripts[name]).read_bytes()).hexdigest()
        statuses = data.setdefault("statuses", {})
        statuses.update({key: "passed" for key in data["results"]})
        statuses.update({key: "capped" for key in data.get("capped", {})})
        for key in KEYS:
            statuses.setdefault(key, "not_implemented" if key not in supported(name) else "not_run")
        path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def supported(name):
    if name in ("sqa", "quantumalgebra"):
        return set(KEYS)
    if name in ("sympy", "sneg"):
        return set(KEYS) - {"tavis_cummings_build", "tavis_cummings_comm", "meanfield_average"}
    return {key for key in KEYS if key.startswith(("fock_", "bh_", "bose_"))}


def validate_dataset(directory):
    run = json.loads((directory / "run.json").read_text())
    expected_jobs = {"julia" if name in ("sqa", "quantumalgebra") else name for name in run["packages"]}
    if run["contract"] != CONTRACT:
        raise ValueError("Wrong run contract")
    if set(run["jobs"]) != expected_jobs or any(value != "passed" for value in run["jobs"].values()):
        raise ValueError("Run has failed or timed-out jobs; inspect its logs")
    for name in run["packages"]:
        data = json.loads((directory / f"{name}.json").read_text())
        if data["meta"]["run_id"] != run["run_id"] or data["meta"]["contract"] != CONTRACT:
            raise ValueError(f"Mixed run/contract: {name}")
        if data["meta"]["contract_sha256"] != run["contract_sha256"]:
            raise ValueError(f"Mixed contract hashes: {name}")
        script = {"sqa": "julia_bench.jl", "quantumalgebra": "julia_bench.jl", "sympy": "sympy_bench.py", "openfermion": "openfermion_bench.py", "sneg": "sneg_bench.wls"}[name]
        if data["meta"]["benchmark_script_sha256"] != hashlib.sha256((directory / "source" / script).read_bytes()).hexdigest():
            raise ValueError(f"Script hash mismatch: {name}")
        expected = supported(name) & set(run["keys"])
        if set(data["results"]) | set(data.get("capped", {})) != expected:
            raise ValueError(f"Missing or unexpected scenarios: {name}")
        if any(not isinstance(value, (int, float)) or not 30 < value < float("inf") for value in data.get("capped", {}).values()):
            raise ValueError(f"Invalid cap measurement: {name}")
        for key, result in data["results"].items():
            if not isinstance(result["time_ns"], (int, float)) or not 0 < result["time_ns"] < float("inf"):
                raise ValueError(f"Invalid timing: {name}/{key}")
    return run


def report(args):
    directory = args.run.resolve()
    run = validate_dataset(directory)
    if args.publish and (set(run["keys"]) != set(KEYS) or not {"sqa", "quantumalgebra", "sympy", "openfermion"}.issubset(run["packages"])):
        raise ValueError("Publication requires the complete core suite, not a smoke run")
    environment = dict(os.environ, SQA_COMPARE_OUTPUT=str(directory),
                       SQA_COMPARE_REPORT=str(directory / "report"),
                       SQA_COMPARE_PUBLISH=str(args.publish).lower())
    subprocess.run([args.julia, "--startup-file=no", "--project=benchmark", str(HERE / "make_table.jl")], cwd=ROOT, env=environment, check=True)
    bundle = directory / "report/reproduction.tar.gz"
    with tarfile.open(bundle, "w:gz") as archive:
        for path in sorted(directory.iterdir()):
            if path.name != "report":
                archive.add(path, arcname=path.name)
        demonstrations = HERE / "demonstrations"
        if demonstrations.exists():
            for source in demonstrations.iterdir():
                if source.suffix in (".jl", ".md", ".toml"):
                    archive.add(source, arcname=f"demonstrations/{source.name}")
            for name in ("high_order", "interoperability"):
                result = HERE / "results" / f"{name}.json"
                if not result.exists():
                    continue
                data = json.loads(result.read_text())
                valid = data.get("verified", False) if name == "interoperability" else bool(data["results"]) and all(row["verified"] for row in data["results"])
                script = demonstrations / f"{name}.jl"
                if not valid or data["meta"]["script_sha256"] != hashlib.sha256(script.read_bytes()).hexdigest():
                    raise ValueError(f"Invalid demonstration artifact: {name}")
                for extension in ("json", "png", "svg"):
                    artifact = result.with_suffix(f".{extension}")
                    if artifact.exists():
                        archive.add(artifact, arcname=f"demonstrations/results/{artifact.name}")
    if args.publish:
        shutil.copy(bundle, HERE / "results/reproduction.tar.gz")
        shutil.copy(directory / "run.json", HERE / "results/run.json")
    print(f"Reproduction archive: {bundle}")


def benchmark(args):
    keys = args.keys.split(",") if args.keys else KEYS
    if set(keys) - set(KEYS):
        raise ValueError("Unknown scenario key")
    run_id = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    directory = (HERE / "runs" / run_id).resolve()
    directory.mkdir(parents=True)
    packages = ["sqa", "quantumalgebra", "sympy", "openfermion"]
    if args.sneg:
        packages.append("sneg")
    context = provenance(__file__)
    context["orchestrator_sha256"] = context.pop("benchmark_script_sha256")
    context.update(run_id=run_id, contract=CONTRACT,
                   contract_sha256=hashlib.sha256((HERE / "BENCHMARKS.md").read_bytes()).hexdigest(),
                   seconds=args.seconds, threads="1", keys=keys)
    source = directory / "source"
    source.mkdir()
    for path in HERE.iterdir():
        if path.suffix in (".jl", ".py", ".wls", ".md", ".in", ".lock"):
            shutil.copy(path, source / path.name)
    for name in ("Project.toml", "Manifest.toml"):
        path = ROOT / "benchmark" / name
        if path.exists():
            shutil.copy(path, directory / name)
    python = args.python or str(HERE / ".venv/bin/python")
    context["julia_runtime"] = subprocess.check_output([args.julia, "--version"], text=True).strip()
    context["python_runtime"] = subprocess.check_output([python, "-c", "import platform; print(platform.python_version())"], text=True).strip()
    frozen = subprocess.check_output(["uv", "pip", "freeze", "--python", python], text=True)
    (directory / "python-requirements.txt").write_text(frozen)
    patch = subprocess.check_output(["git", "-C", str(ROOT), "diff", "HEAD", "--", "src", "ext", "Project.toml"], text=True)
    (directory / "source.patch").write_text(patch)
    manifest = dict(context, packages=packages, jobs={})
    path = directory / "run.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    with tempfile.TemporaryDirectory(prefix="sqa-compare-") as temporary:
        environment = dict(os.environ, SQA_COMPARE_ROOT=str(ROOT), SQA_COMPARE_OUTPUT=str(directory), SQA_COMPARE_TEMP=temporary,
                           SQA_COMPARE_RUN_ID=run_id, SQA_COMPARE_KEYS=",".join(keys),
                           SQA_COMPARE_SECONDS=str(args.seconds), JULIA_NUM_THREADS="1",
                           OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        environment["SQA_COMPARE_VALIDATE_ONLY"] = str(args.validate_only).lower()
        jobs = [("julia", [args.julia, "--startup-file=no", "--threads=1", "--gcthreads=1", "--project=benchmark", str(source / "julia_bench.jl")]),
                ("sympy", [python, str(source / "sympy_bench.py")]),
                ("openfermion", [python, str(source / "openfermion_bench.py")])]
        if args.sneg:
            sneg_source = Path(environment.get("SNEG_PATH", "")) / "sneg.m"
            if sneg_source.is_file():
                shutil.copy(sneg_source, source / "sneg.m")
                environment["SNEG_PATH"] = str(source)
                manifest["sneg_source_sha256"] = hashlib.sha256(sneg_source.read_bytes()).hexdigest()
            jobs.append(("sneg", ["wolframscript", "-f", str(source / "sneg_bench.wls")]))
        for name, command in jobs:
            print(f"{name}: running; log {directory / (name + '.log')}", flush=True)
            try:
                status = execute(command, environment, directory / f"{name}.log", args.timeout)
            except FileNotFoundError as error:
                (directory / f"{name}.log").write_text(str(error))
                status = "not_installed"
            manifest["jobs"][name] = status
            path.write_text(json.dumps(manifest, indent=2) + "\n")
            print(f"{name}: {status}", flush=True)
    if any(value != "passed" for value in manifest["jobs"].values()):
        raise RuntimeError(f"Run incomplete: {directory}; see run.json and logs")
    if not args.validate_only:
        enrich(directory, packages, context)
        validate_dataset(directory)
    print(directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia", default="julia")
    commands = parser.add_subparsers(dest="command", required=True)
    setup = commands.add_parser("prepare")
    modes = setup.add_mutually_exclusive_group()
    modes.add_argument("--refresh", action="store_true", help="Intentionally update dependency locks")
    modes.add_argument("--from-run", type=Path, help="Restore resolved environments from a run archive")
    run = commands.add_parser("run")
    run.add_argument("--python", help="Python interpreter with locked dependencies")
    run.add_argument("--keys", default="", help="Comma-separated subset for smoke runs")
    run.add_argument("--seconds", type=float, default=10)
    run.add_argument("--timeout", type=float, default=1800, help="Hard wall deadline for each package job, including warmup")
    run.add_argument("--sneg", action="store_true", help="Include optional locally installed Mathematica/SNEG")
    run.add_argument("--validate-only", action="store_true")
    render = commands.add_parser("report")
    render.add_argument("run", type=Path)
    render.add_argument("--publish", action="store_true")
    args = parser.parse_args()
    if getattr(args, "seconds", 1) <= 0 or getattr(args, "timeout", 1) <= 0:
        parser.error("Timing budgets must be positive")
    {"prepare": prepare, "run": benchmark, "report": report}[args.command](args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Negative controls: planted false identities must be detected, with and without -O.

For each shipped checker, a copy with exactly one falsified expected value is
written to a fresh working directory and run as ``python -I`` and as
``python -I -O``. A control passes only if the mutant exits nonzero with the
explicit "CHECK FAILED" message (not, for example, a syntax error). Unmutated
fast checkers are also run under ``-I -O`` as positive controls.

Standard library only; it runs the copies with the current interpreter, so
NumPy/SciPy/mpmath must be importable for the numerical scripts.
Writes negative_controls.json to --out DIR (default: current directory).
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

# (script, original text, falsified text, extra arguments, description)
MUTATIONS = [
    ("verify_gw.py", "== 24*(ai*ai+aj*aj)-32*ai*aj", "== 25*(ai*ai+aj*aj)-32*ai*aj",
     ["--exact-only"], "441-case grid: generator norm 24 -> 25"),
    ("independent_symbolic_audit.py", "24 * (ai**2 + aj**2) - 32 * ai * aj",
     "25 * (ai**2 + aj**2) - 32 * ai * aj", [], "generator norm identity 24 -> 25"),
    ("verify_combined_flow.py", "check(K <= 32*A*B*S+1e-7", "check(K <= 3*A*B*S+1e-7",
     [], "curvature bound K_C <= 32 A2 B2 S -> 3 A2 B2 S"),
    ("proof_digest_checks.py", "for c in deterministic)==8,", "for c in deterministic)==9,",
     [], "finite-square optimum 8 -> 9"),
    ("extension_checks.py", "check(QG == 4 * (sum(a) - sum(b)) ** 2 + 8 * nab",
     "check(QG == 4 * (sum(a) - sum(b)) ** 2 + 9 * nab", [], "Q_G identity coefficient 8 -> 9"),
    ("gw_gaussian_bracket.py", '"sharp_32S2_over_KC": 32 * S * S / k_c',
     '"sharp_32S2_over_KC": 33 * S * S / k_c', [], "sharp certificate 32 -> 33 (demo self-check)"),
]
SUPPORT = {"extension_checks.py": ["gw_gaussian_bracket.py"]}
POSITIVE = [("independent_symbolic_audit.py", []), ("verify_gw.py", ["--exact-only"]),
            ("gw_gaussian_bracket.py", [])]
FLAGS = {"-I": ["-I"], "-I -O": ["-I", "-O"]}


def check(condition, message):
    if not condition:
        raise SystemExit(f"CHECK FAILED: {message}")


def run(script, args, flags, out_dir, timeout=900):
    command = [sys.executable, *flags, str(script), *args]
    if script.name in ("verify_gw.py", "verify_combined_flow.py", "proof_digest_checks.py",
                       "extension_checks.py"):
        command += ["--out", str(out_dir)]
    started = time.perf_counter()
    proc = subprocess.run(command, capture_output=True, text=True, timeout=timeout, cwd=out_dir)
    return proc, time.perf_counter() - started


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=".", help="directory for negative_controls.json")
    parser.add_argument("--work", default=None, help="new working directory (default: fresh temp dir)")
    args = parser.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    work = Path(args.work) if args.work else Path(tempfile.mkdtemp(prefix="gw-negative-controls-"))
    work.mkdir(parents=True, exist_ok=True)
    check(not any(work.iterdir()), f"working directory {work} must be empty")
    records = []
    for index, (name, old, new, extra, description) in enumerate(MUTATIONS):
        source = (HERE / name).read_text()
        check(source.count(old) == 1, f"mutation anchor not unique in {name}: {old!r}")
        case_dir = work / f"mutant_{index}_{Path(name).stem}"
        case_dir.mkdir()
        (case_dir / name).write_text(source.replace(old, new))
        for support in SUPPORT.get(name, []):
            (case_dir / support).write_text((HERE / support).read_text())
        for label, flags in FLAGS.items():
            run_dir = case_dir / ("run_O" if "-O" in flags else "run")
            run_dir.mkdir()
            proc, seconds = run(case_dir / name, extra, flags, run_dir)
            detected = proc.returncode != 0 and "CHECK FAILED" in proc.stderr
            message = next((line for line in proc.stderr.splitlines() if "CHECK FAILED" in line), "")
            records.append({"script": name, "mutation": description, "flags": label,
                            "exit_code": proc.returncode, "detected": detected,
                            "message": message[:300], "seconds": round(seconds, 2)})
            check(detected, f"mutant of {name} not detected under python {label}: "
                            f"exit {proc.returncode}, stderr tail {proc.stderr[-300:]!r}")
    positives = []
    for name, extra in POSITIVE:
        run_dir = work / f"positive_{Path(name).stem}"
        run_dir.mkdir()
        proc, seconds = run(HERE / name, extra, FLAGS["-I -O"], run_dir)
        positives.append({"script": name, "flags": "-I -O", "exit_code": proc.returncode,
                          "seconds": round(seconds, 2)})
        check(proc.returncode == 0, f"unmutated {name} failed under -I -O: {proc.stderr[-300:]!r}")
    result = {"python": platform.python_version(), "interpreter": sys.executable,
              "working_directory": str(work), "negative_controls": records,
              "positive_controls": positives,
              "summary": f"{sum(r['detected'] for r in records)}/{len(records)} planted false "
                         "identities detected (each under python -I and python -I -O)"}
    # Status is assigned only after every control has been checked.
    result["status"] = "PASS"
    (out / "negative_controls.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "summary": result["summary"]}, indent=2))


if __name__ == "__main__":
    main()

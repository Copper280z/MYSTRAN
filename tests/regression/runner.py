#!/usr/bin/env python3
"""
Regression test runner for MYSTRAN.

Runs two MYSTRAN binaries (baseline and current) on a curated set of benchmark
cases and compares their F06 output files using f06diff.  Any numerical
differences exceeding the thresholds cause a test failure.

Usage:
    python runner.py \\
        --baseline  <path/to/mystran_baseline> \\
        --current   <path/to/mystran_current> \\
        --f06diff   <path/to/f06diff> \\
        --decks-dir <path/to/MYSTRAN_Benchmark/Benchmark_Decks> \\
        [--cases    <path/to/benchmark_cases.txt>] \\
        [--workdir  <scratch_dir>] \\
        [--difference 1e-6] \\
        [--ratio     1e-4]

Exit code: 0 if all cases match within thresholds, 1 otherwise.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_CASES_FILE = Path(__file__).parent / "benchmark_cases.txt"
DEFAULT_DIFFERENCE = 1e-6
DEFAULT_RATIO = 1e-4


def run_mystran(mystran_exe: Path, dat_file: Path, run_dir: Path) -> bool:
    """Copy dat_file into run_dir and run mystran on it.  Returns True on success."""
    dest = run_dir / dat_file.name
    shutil.copy(dat_file, dest)
    result = subprocess.run(
        [str(mystran_exe), dest.name],
        cwd=str(run_dir),
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


def find_f06(run_dir: Path, stem: str) -> Path | None:
    """Locate the F06 output file produced by MYSTRAN (may be .F06 or .f06)."""
    for ext in (".F06", ".f06"):
        p = run_dir / (stem + ext)
        if p.exists():
            return p
    return None


def compare_f06(f06diff_exe: Path, baseline_f06: Path, current_f06: Path,
                difference: float, ratio: float) -> tuple[bool, str]:
    """
    Run f06diff and return (passed, output_text).
    f06diff exits 0 when files match within thresholds, non-zero otherwise.
    """
    cmd = [
        str(f06diff_exe),
        "--difference", str(difference),
        "--ratio",      str(ratio),
        str(baseline_f06),
        str(current_f06),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = (result.stdout + result.stderr).strip()
    return result.returncode == 0, output


def load_cases(cases_file: Path) -> list[str]:
    cases = []
    with open(cases_file) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                cases.append(line)
    return cases


def main():
    parser = argparse.ArgumentParser(description="MYSTRAN regression test runner")
    parser.add_argument("--baseline",   required=True, help="Baseline MYSTRAN binary")
    parser.add_argument("--current",    required=True, help="Current MYSTRAN binary to test")
    parser.add_argument("--f06diff",    required=True, help="f06diff binary")
    parser.add_argument("--decks-dir",  required=True, help="Directory containing benchmark .DAT files")
    parser.add_argument("--cases",      default=str(DEFAULT_CASES_FILE),
                        help=f"Text file listing cases to run (default: {DEFAULT_CASES_FILE})")
    parser.add_argument("--workdir",    default=None,
                        help="Scratch directory (created if absent; default: temp dir)")
    parser.add_argument("--difference", type=float, default=DEFAULT_DIFFERENCE,
                        help=f"Max absolute difference threshold (default: {DEFAULT_DIFFERENCE})")
    parser.add_argument("--ratio",      type=float, default=DEFAULT_RATIO,
                        help=f"Max big-to-small ratio threshold (default: {DEFAULT_RATIO})")
    args = parser.parse_args()

    baseline_exe  = Path(args.baseline).resolve()
    current_exe   = Path(args.current).resolve()
    f06diff_exe   = Path(args.f06diff).resolve()
    decks_dir     = Path(args.decks_dir).resolve()
    cases_file    = Path(args.cases).resolve()

    for p, label in [(baseline_exe, "baseline"), (current_exe, "current"),
                     (f06diff_exe, "f06diff"), (decks_dir, "decks-dir"),
                     (cases_file, "cases file")]:
        if not p.exists():
            sys.exit(f"ERROR: {label} not found: {p}")

    cases = load_cases(cases_file)
    if not cases:
        sys.exit(f"ERROR: no cases found in {cases_file}")

    use_temp = args.workdir is None
    workdir = Path(tempfile.mkdtemp(prefix="mystran_regression_")) if use_temp \
              else Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    print(f"Baseline:  {baseline_exe}")
    print(f"Current:   {current_exe}")
    print(f"f06diff:   {f06diff_exe}")
    print(f"Workdir:   {workdir}")
    print(f"Cases:     {len(cases)}")
    print(f"Threshold: difference={args.difference:.1e}  ratio={args.ratio:.1e}")
    print()

    n_pass = n_fail = n_error = 0

    for case_rel in cases:
        dat_file = decks_dir / case_rel
        if not dat_file.exists():
            print(f"ERROR  {case_rel}: .DAT file not found at {dat_file}")
            n_error += 1
            continue

        stem = dat_file.stem
        base_run = workdir / "baseline" / stem
        curr_run  = workdir / "current"  / stem
        base_run.mkdir(parents=True, exist_ok=True)
        curr_run.mkdir(parents=True, exist_ok=True)

        # Run both binaries
        ok_base = run_mystran(baseline_exe, dat_file, base_run)
        ok_curr = run_mystran(current_exe,  dat_file, curr_run)

        if not ok_base:
            print(f"ERROR  {case_rel}: baseline run failed")
            n_error += 1
            continue
        if not ok_curr:
            print(f"ERROR  {case_rel}: current run failed")
            n_error += 1
            continue

        # Locate F06 outputs
        f06_base = find_f06(base_run, stem)
        f06_curr = find_f06(curr_run,  stem)

        if f06_base is None:
            print(f"ERROR  {case_rel}: baseline F06 not found")
            n_error += 1
            continue
        if f06_curr is None:
            print(f"ERROR  {case_rel}: current F06 not found")
            n_error += 1
            continue

        passed, diff_output = compare_f06(f06diff_exe, f06_base, f06_curr,
                                          args.difference, args.ratio)
        if passed:
            print(f"PASS   {case_rel}")
            n_pass += 1
        else:
            print(f"FAIL   {case_rel}")
            for line in diff_output.splitlines()[:20]:   # cap output
                print(f"  {line}")
            if diff_output.count("\n") > 20:
                print(f"  ... (truncated)")
            n_fail += 1

    if use_temp and n_fail == 0 and n_error == 0:
        shutil.rmtree(workdir, ignore_errors=True)

    print()
    print(f"Total: {n_pass} passed, {n_fail} failed, {n_error} error(s)")

    if n_fail > 0 or n_error > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

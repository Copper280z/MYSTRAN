#!/usr/bin/env python3
"""
E2E test runner for MYSTRAN.

For each *.bdf file in the given cases directory, runs MYSTRAN, parses the
resulting OP2 file with pyNastran, and compares values against the matching
*.expected.json file.

Usage:
    python runner.py <mystran_binary> <cases_dir>

Exit code: 0 if all checks pass, 1 if any fail or any case errors out.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

try:
    from pyNastran.op2.op2 import OP2
except ImportError:
    sys.exit("ERROR: pyNastran is not installed. Run: pip install pyNastran")

# ---------------------------------------------------------------------------
# DOF name → 0-based column index in pyNastran displacement data
# ---------------------------------------------------------------------------
DOF_INDEX = {"T1": 0, "T2": 1, "T3": 2, "R1": 3, "R2": 4, "R3": 5}

# ---------------------------------------------------------------------------
# Shell stress component → 0-based column index (MSC Nastran convention)
# Columns: fiber_dist, sx, sy, txy, angle, omax, omin, von_mises
# ---------------------------------------------------------------------------
SHELL_STRESS_COL = {
    "fd": 0,          # fiber distance
    "sx": 1,
    "sy": 2,
    "sxy": 3,
    "txy": 3,
    "angle": 4,
    "omax": 5,
    "omin": 6,
    "von_mises": 7,
    "max_shear": 7,
}

# Solid element stress components — layout varies by element; pyNastran stores
# [sx, sy, sz, txy, tyz, txz, omax, omid, omin, von_mises]  for CHEXA/CTETRA
SOLID_STRESS_COL = {
    "sx": 0,
    "sy": 1,
    "sz": 2,
    "sxy": 3, "txy": 3,
    "syz": 4, "tyz": 4,
    "sxz": 5, "txz": 5,
    "omax": 6,
    "omid": 7,
    "omin": 8,
    "von_mises": 9,
}

# Table attribute names to search for stress results in the OP2 object.
# Ordered from most to least specific so the first match wins.
STRESS_TABLE_NAMES = [
    "cquad4_stress", "ctria3_stress", "ctria6_stress", "cquad8_stress",
    "chexa_stress", "cpenta_stress", "ctetra_stress",
    "cbar_stress", "cbeam_stress", "crod_stress",
]


def _compare(actual, expected, rel_tol=None, abs_tol=None, label=""):
    """Return (passed, message).  Exactly one of rel_tol / abs_tol must be set."""
    if rel_tol is not None:
        if expected == 0.0:
            passed = abs(actual) <= rel_tol
        else:
            passed = abs(actual - expected) / abs(expected) <= rel_tol
        msg = (f"{label}: actual={actual:.6g}  expected={expected:.6g}  "
               f"rel_err={abs(actual-expected)/max(abs(expected),1e-300):.3e}  "
               f"rel_tol={rel_tol:.3e}")
    else:
        passed = abs(actual - expected) <= abs_tol
        msg = (f"{label}: actual={actual:.6g}  expected={expected:.6g}  "
               f"abs_err={abs(actual-expected):.3e}  abs_tol={abs_tol:.3e}")
    return passed, msg


def _find_node_index(obj, node_id):
    """Return the row index for node_id in a displacement/eigenvector result."""
    import numpy as np
    # node_gridtype shape: (nnodes, 2); column 0 is node ID
    matches = np.where(obj.node_gridtype[:, 0] == node_id)[0]
    if len(matches) == 0:
        return None
    return int(matches[0])


def _find_element_row(obj, elem_id):
    """
    Return the row index for elem_id in a stress result object.
    Shell elements: obj.element_node[:, 0] contains element IDs (possibly
    repeated for each fibre/layer).  Solid elements: obj.element[:] .
    """
    import numpy as np
    if hasattr(obj, "element_node"):
        # Shell: (nrows, 2) — column 0 is elem ID; take the first centroid row
        matches = np.where(obj.element_node[:, 0] == elem_id)[0]
        if len(matches) > 0:
            return int(matches[0])
    if hasattr(obj, "element"):
        matches = np.where(obj.element == elem_id)[0]
        if len(matches) > 0:
            return int(matches[0])
    return None


def _is_shell_stress(table_name):
    return any(k in table_name for k in ("quad", "tria", "bar", "beam", "rod"))


def check_node_displacement(op2, subcase, node_id, dof, expected, rel_tol, abs_tol):
    """Look up a nodal displacement and compare."""
    disp_tables = {}
    for attr in ("displacements", "eigenvectors"):
        tbl = getattr(op2, attr, None)
        if tbl:
            disp_tables.update(tbl)

    if subcase not in disp_tables:
        return False, f"  subcase {subcase} not found in displacement/eigenvector results"

    obj = disp_tables[subcase]
    inode = _find_node_index(obj, node_id)
    if inode is None:
        return False, f"  node {node_id} not found in result"

    idof = DOF_INDEX.get(dof)
    if idof is None:
        return False, f"  unknown DOF '{dof}'"

    actual = float(obj.data[0, inode, idof])
    return _compare(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol,
                    label=f"node {node_id} {dof}")


def check_element_stress(op2, subcase, elem_id, component, expected, rel_tol, abs_tol):
    """Search all stress tables for elem_id and compare the requested component."""
    for table_name in STRESS_TABLE_NAMES:
        tbl = getattr(op2, table_name, None)
        if not tbl or subcase not in tbl:
            continue
        obj = tbl[subcase]
        irow = _find_element_row(obj, elem_id)
        if irow is None:
            continue

        # Pick the right column map
        if _is_shell_stress(table_name):
            col_map = SHELL_STRESS_COL
        else:
            col_map = SOLID_STRESS_COL

        icol = col_map.get(component)
        if icol is None:
            return False, f"  unknown stress component '{component}' for {table_name}"

        actual = float(obj.data[0, irow, icol])
        return _compare(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol,
                        label=f"elem {elem_id} {component}")

    return False, f"  elem {elem_id} not found in any stress table for subcase {subcase}"


def check_eigenfrequency(op2, subcase, mode, expected, rel_tol, abs_tol):
    """Look up an eigenfrequency (Hz) for the given mode number (1-indexed)."""
    # pyNastran stores eigenvalues keyed by (subcase_id, 'LAMA') or plain subcase_id
    eig_tables = getattr(op2, "eigenvalues", None) or {}
    obj = None
    for key, val in eig_tables.items():
        sc = key[0] if isinstance(key, tuple) else key
        if sc == subcase:
            obj = val
            break

    if obj is None:
        # Fall back: check eigenvectors table for frequency data
        evec = (getattr(op2, "eigenvectors", None) or {}).get(subcase)
        if evec is not None and hasattr(evec, "freq"):
            import numpy as np
            freqs = evec.freq
            if mode - 1 < len(freqs):
                actual = float(freqs[mode - 1])
                return _compare(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol,
                                label=f"mode {mode} freq (Hz)")
        return False, f"  eigenvalue table not found for subcase {subcase}"

    import numpy as np
    if hasattr(obj, "frequency"):
        freqs = obj.frequency
    elif hasattr(obj, "Hz"):
        freqs = obj.Hz
    else:
        return False, f"  eigenvalue object has no frequency attribute"

    if mode - 1 >= len(freqs):
        return False, f"  mode {mode} out of range (only {len(freqs)} modes)"

    actual = float(freqs[mode - 1])
    return _compare(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol,
                    label=f"mode {mode} freq (Hz)")


def run_case(mystran_exe, bdf_path, expected_json_path):
    """
    Run MYSTRAN on bdf_path, parse the OP2, check all entries in expected_json.
    Returns (n_passed, n_failed, n_errors, messages).
    """
    bdf_path = Path(bdf_path).resolve()
    op2_path = bdf_path.with_suffix(".OP2")
    if not op2_path.exists():
        op2_path = bdf_path.with_suffix(".op2")

    # Run MYSTRAN
    result = subprocess.run(
        [str(mystran_exe), str(bdf_path.name)],
        cwd=str(bdf_path.parent),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return 0, 0, 1, [f"  MYSTRAN exited with code {result.returncode}"]

    # Locate OP2
    op2_path = bdf_path.with_suffix(".OP2")
    if not op2_path.exists():
        op2_path = bdf_path.with_suffix(".op2")
    if not op2_path.exists():
        return 0, 0, 1, [f"  OP2 file not found: {bdf_path.stem}.op2"]

    # Parse OP2
    try:
        op2 = OP2(debug=False)
        op2.read_op2(str(op2_path))
    except Exception as exc:
        return 0, 0, 1, [f"  Failed to parse OP2: {exc}"]

    # Load expected results
    with open(expected_json_path) as fh:
        spec = json.load(fh)

    n_pass = 0
    n_fail = 0
    messages = []

    for chk in spec.get("checks", []):
        chk_type = chk["type"]
        subcase = chk.get("subcase", 1)
        expected = chk["expected"]
        rel_tol = chk.get("rel_tol", None)
        abs_tol = chk.get("abs_tol", None)
        if rel_tol is None and abs_tol is None:
            abs_tol = 1e-8  # fallback

        try:
            if chk_type == "node_displacement":
                passed, msg = check_node_displacement(
                    op2, subcase, chk["node_id"], chk["dof"],
                    expected, rel_tol, abs_tol)
            elif chk_type == "element_stress":
                passed, msg = check_element_stress(
                    op2, subcase, chk["elem_id"], chk["component"],
                    expected, rel_tol, abs_tol)
            elif chk_type == "eigenfrequency":
                passed, msg = check_eigenfrequency(
                    op2, subcase, chk["mode"],
                    expected, rel_tol, abs_tol)
            else:
                passed, msg = False, f"  unknown check type '{chk_type}'"
        except Exception as exc:
            passed, msg = False, f"  exception during check: {exc}"

        if passed:
            n_pass += 1
        else:
            n_fail += 1
            comment = chk.get("comment", "")
            messages.append(f"  FAIL [{chk_type}] {comment}")
            messages.append(f"    {msg}")

    return n_pass, n_fail, 0, messages


def main():
    parser = argparse.ArgumentParser(description="MYSTRAN E2E test runner")
    parser.add_argument("mystran", help="Path to the mystran binary")
    parser.add_argument("cases_dir", help="Directory containing .bdf and .expected.json files")
    args = parser.parse_args()

    mystran_exe = Path(args.mystran).resolve()
    if not mystran_exe.exists():
        sys.exit(f"ERROR: mystran binary not found: {mystran_exe}")

    cases_dir = Path(args.cases_dir).resolve()
    bdf_files = sorted(cases_dir.glob("*.bdf"))
    if not bdf_files:
        sys.exit(f"ERROR: no .bdf files found in {cases_dir}")

    total_pass = total_fail = total_error = 0

    for bdf_path in bdf_files:
        json_path = bdf_path.with_suffix(".expected.json")
        if not json_path.exists():
            print(f"SKIP {bdf_path.name}: no matching .expected.json")
            continue

        n_pass, n_fail, n_err, messages = run_case(mystran_exe, bdf_path, json_path)
        total_pass += n_pass
        total_fail += n_fail
        total_error += n_err

        status = "PASS" if (n_fail == 0 and n_err == 0) else "FAIL"
        summary = f"{n_pass} passed"
        if n_fail:
            summary += f", {n_fail} failed"
        if n_err:
            summary += f", {n_err} error(s)"
        print(f"{status} {bdf_path.name}  [{summary}]")
        for msg in messages:
            print(msg)

    print()
    print(f"Total: {total_pass} passed, {total_fail} failed, {total_error} error(s)")

    if total_fail > 0 or total_error > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

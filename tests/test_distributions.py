#!/usr/bin/env python3
"""
Statistical distribution tests for FARM simulation.
Generates events and compares physics distributions against reference
using Kolmogorov-Smirnov tests.

Usage:
  # Generate reference:
  python3 test_distributions.py --sim ./builddir/simulation --workdir . --generate-ref

  # Run comparison test:
  python3 test_distributions.py --sim ./builddir/simulation --workdir .
"""

import sys
import os
import argparse
import subprocess
import tempfile
import numpy as np

# Add scripts directory to path for importing plot_lund
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts'))
from plot_lund import read_lund, compute_dis


REF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "statistical", "ref")
NEVENT = 50000
SEED = 123
TARGET = 1  # deuterium
NKIN = 2000
E0 = 10.5
ALPHA = 0.001  # p-value threshold


def run_simulation(sim_binary, output_file, workdir):
    """Run the FARM simulation with fixed parameters."""
    cmd = [
        sim_binary,
        "--nevent", str(NEVENT),
        "--output", output_file,
        "--target", str(TARGET),
        "--nkin", str(NKIN),
        "--e0", str(E0),
        "--seed", str(SEED),
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    if result.returncode != 0:
        print("STDERR:", result.stderr)
        raise RuntimeError(f"Simulation failed with exit code {result.returncode}")
    print(f"Simulation completed: {NEVENT} events")


def extract_distributions(lund_file):
    """Extract physics distributions from a LUND file."""
    events = read_lund(lund_file)
    dis = compute_dis(events)

    # Hadron distributions
    z_vals = []
    pt2_vals = []
    multiplicities = []

    for ev in events:
        nhad = 0
        for p in ev["particles"]:
            pid = abs(p["pid"])
            # Count charged hadrons (pions, kaons, protons)
            if pid in (211, 321, 2212):
                nhad += 1
                px, py, pz, E = p["px"], p["py"], p["pz"], p["E"]
                p_mag = np.sqrt(px**2 + py**2 + pz**2)
                pt2 = px**2 + py**2
                pt2_vals.append(pt2)

                # Find scattered electron for z calculation
                e_beam = ev["beam_energy"]
                for ep in ev["particles"]:
                    if ep["pid"] == 11:
                        nu = e_beam - ep["E"]
                        if nu > 0:
                            z = E / nu
                            z_vals.append(z)
                        break
        multiplicities.append(nhad)

    return {
        "Q2": np.array(dis["Q2"]),
        "xB": np.array(dis["xB"]),
        "W": np.array(dis["W"]),
        "nu": np.array(dis["nu"]),
        "y": np.array(dis["y"]),
        "z": np.array(z_vals),
        "pt2": np.array(pt2_vals),
        "multiplicity": np.array(multiplicities),
    }


def generate_reference(sim_binary, workdir):
    """Generate reference distributions and save as .npz."""
    os.makedirs(REF_DIR, exist_ok=True)
    ref_file = os.path.join(REF_DIR, "dis_distributions.npz")

    with tempfile.NamedTemporaryFile(suffix=".lund", delete=False) as f:
        tmp_lund = f.name

    try:
        run_simulation(sim_binary, tmp_lund, workdir)
        dists = extract_distributions(tmp_lund)
        np.savez(ref_file, **dists)
        print(f"Reference saved to {ref_file}")
        for name, arr in dists.items():
            print(f"  {name}: {len(arr)} entries, mean={np.mean(arr):.4f}")
    finally:
        os.unlink(tmp_lund)


def run_ks_test(ref, test, name):
    """Run a 2-sample KS test and return (statistic, p-value)."""
    from scipy.stats import ks_2samp
    stat, pvalue = ks_2samp(ref, test)
    return stat, pvalue


def run_comparison(sim_binary, workdir):
    """Run simulation and compare distributions against reference."""
    ref_file = os.path.join(REF_DIR, "dis_distributions.npz")
    if not os.path.exists(ref_file):
        print(f"Reference file not found: {ref_file}")
        print("Generate it first with --generate-ref")
        sys.exit(2)

    ref_data = np.load(ref_file)

    with tempfile.NamedTemporaryFile(suffix=".lund", delete=False) as f:
        tmp_lund = f.name

    try:
        run_simulation(sim_binary, tmp_lund, workdir)
        test_data = extract_distributions(tmp_lund)
    finally:
        os.unlink(tmp_lund)

    # Run KS tests
    failures = 0
    variables = ["Q2", "xB", "W", "nu", "y", "z", "pt2", "multiplicity"]

    print("\n=== Kolmogorov-Smirnov Test Results ===")
    print(f"{'Variable':<15} {'KS stat':<12} {'p-value':<12} {'Result'}")
    print("-" * 55)

    for var in variables:
        if var not in ref_data or var not in test_data:
            print(f"{var:<15} {'SKIP':>36} (missing data)")
            continue

        ref_arr = ref_data[var]
        test_arr = test_data[var]

        if len(ref_arr) == 0 or len(test_arr) == 0:
            print(f"{var:<15} {'SKIP':>36} (empty array)")
            continue

        stat, pvalue = run_ks_test(ref_arr, test_arr, var)
        passed = pvalue > ALPHA
        status = "PASS" if passed else "FAIL"

        print(f"{var:<15} {stat:<12.6f} {pvalue:<12.6f} {status}")

        if not passed:
            failures += 1

    print("-" * 55)
    if failures == 0:
        print(f"All distribution tests PASSED (alpha={ALPHA})")
        return 0
    else:
        print(f"{failures} distribution test(s) FAILED")
        return 1


def main():
    parser = argparse.ArgumentParser(description="Statistical distribution tests for FARM")
    parser.add_argument("--sim", required=True, help="Path to simulation binary")
    parser.add_argument("--workdir", default=".", help="Working directory for simulation")
    parser.add_argument("--generate-ref", action="store_true",
                        help="Generate reference distributions instead of testing")
    args = parser.parse_args()

    if args.generate_ref:
        generate_reference(args.sim, args.workdir)
        sys.exit(0)
    else:
        sys.exit(run_comparison(args.sim, args.workdir))


if __name__ == "__main__":
    main()

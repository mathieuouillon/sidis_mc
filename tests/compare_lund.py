#!/usr/bin/env python3
"""
Tolerance-based comparison of two LUND format files.
Used as a fallback when exact diff fails (e.g., cross-platform).

Usage: python3 compare_lund.py <reference.lund> <test.lund> [--rel-tol 1e-6] [--abs-tol 1e-10]
"""

import sys
import argparse


def parse_lund_file(filename):
    """Parse a LUND file into a list of events, each with header and particles."""
    events = []
    with open(filename) as f:
        while True:
            header = f.readline()
            if not header or not header.strip():
                break
            hfields = header.split()
            if len(hfields) < 10:
                break
            npart = int(hfields[0])
            particles = []
            for _ in range(npart):
                pline = f.readline()
                if not pline:
                    break
                particles.append(pline.split())
            events.append({"header": hfields, "particles": particles})
    return events


def floats_match(a_str, b_str, rel_tol, abs_tol):
    """Compare two string representations of numbers with tolerance."""
    try:
        a = float(a_str)
        b = float(b_str)
    except ValueError:
        return a_str == b_str

    if a == b:
        return True
    diff = abs(a - b)
    if diff <= abs_tol:
        return True
    denom = max(abs(a), abs(b))
    if denom > 0 and diff / denom <= rel_tol:
        return True
    return False


def compare_events(ref_events, test_events, rel_tol, abs_tol):
    """Compare two lists of LUND events. Returns (passed, message)."""
    if len(ref_events) != len(test_events):
        return False, f"Event count mismatch: ref={len(ref_events)}, test={len(test_events)}"

    for iev, (ref_ev, test_ev) in enumerate(zip(ref_events, test_events)):
        ref_h = ref_ev["header"]
        test_h = test_ev["header"]

        # Compare header: integer fields exact, float fields with tolerance
        # Header format: Npart(int) iA(int) iZ(int) 0(int) 0(int) 11(int) E0(float) nucleon(int) 1(int) 1.0(float)
        header_int_indices = [0, 1, 2, 3, 4, 5, 7, 8]
        header_float_indices = [6, 9]

        for idx in header_int_indices:
            if idx < len(ref_h) and idx < len(test_h):
                if ref_h[idx] != test_h[idx]:
                    return False, (
                        f"Event {iev+1} header field {idx}: "
                        f"ref={ref_h[idx]}, test={test_h[idx]}"
                    )

        for idx in header_float_indices:
            if idx < len(ref_h) and idx < len(test_h):
                if not floats_match(ref_h[idx], test_h[idx], rel_tol, abs_tol):
                    return False, (
                        f"Event {iev+1} header field {idx}: "
                        f"ref={ref_h[idx]}, test={test_h[idx]}"
                    )

        # Compare particles
        ref_parts = ref_ev["particles"]
        test_parts = test_ev["particles"]
        if len(ref_parts) != len(test_parts):
            return False, (
                f"Event {iev+1} particle count: "
                f"ref={len(ref_parts)}, test={len(test_parts)}"
            )

        for ipart, (ref_p, test_p) in enumerate(zip(ref_parts, test_parts)):
            # Particle format: idx(int) charge(int) 1(int) pid(int) mother(int) 0(int)
            #                  px(float) py(float) pz(float) E(float) mass(float)
            #                  0(int) 0(int) vz(float)
            part_int_indices = [0, 1, 2, 3, 4, 5, 11, 12]
            part_float_indices = [6, 7, 8, 9, 10, 13]

            for idx in part_int_indices:
                if idx < len(ref_p) and idx < len(test_p):
                    if ref_p[idx] != test_p[idx]:
                        return False, (
                            f"Event {iev+1} particle {ipart+1} field {idx}: "
                            f"ref={ref_p[idx]}, test={test_p[idx]}"
                        )

            for idx in part_float_indices:
                if idx < len(ref_p) and idx < len(test_p):
                    if not floats_match(ref_p[idx], test_p[idx], rel_tol, abs_tol):
                        return False, (
                            f"Event {iev+1} particle {ipart+1} field {idx}: "
                            f"ref={ref_p[idx]}, test={test_p[idx]}"
                        )

    return True, f"All {len(ref_events)} events match within tolerance"


def main():
    parser = argparse.ArgumentParser(description="Compare two LUND files with tolerance")
    parser.add_argument("reference", help="Reference LUND file")
    parser.add_argument("test", help="Test LUND file to compare")
    parser.add_argument("--rel-tol", type=float, default=1e-6, help="Relative tolerance")
    parser.add_argument("--abs-tol", type=float, default=1e-10, help="Absolute tolerance")
    args = parser.parse_args()

    ref_events = parse_lund_file(args.reference)
    test_events = parse_lund_file(args.test)

    passed, message = compare_events(ref_events, test_events, args.rel_tol, args.abs_tol)

    if passed:
        print(f"PASS: {message}")
        sys.exit(0)
    else:
        print(f"FAIL: {message}")
        sys.exit(1)


if __name__ == "__main__":
    main()

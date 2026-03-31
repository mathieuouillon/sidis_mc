#!/bin/bash
# Golden file regression test for FARM simulation
# Usage: ./test_golden.sh <config_name> <simulation_binary>
#
# Each config uses --seed 42 for reproducibility and the correct iFM model:
#   iFM=0: proton (rFM=0, no Fermi motion)
#   iFM=5: deuterium, tritium, he3, he4, li6, li7 (R. Wiringa tables)
#   iFM=1: carbon, aluminum, iron, tin, lead (Fermi sphere + Bodek-Ritchie tail)
#   iFM=2: neon, krypton, xenon, copper (Accardi SVG)
#   iFM=3: iron_cs (Accardi CS)
#   iFM=4: carbon_hs (hard sphere)

set -euo pipefail

CONFIG="${1:-proton}"
SIM="${2:-./build/simulation}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GOLDEN_DIR="${SCRIPT_DIR}/golden"
TMP_DIR=$(mktemp -d)
trap "rm -rf ${TMP_DIR}" EXIT

# Define test configurations
case "${CONFIG}" in
    # ── iFM=0: no Fermi motion (rFM=0) ──────────────────────────────────
    proton)
        # target 0: rFM=0, no FM needed
        ARGS="--nevent 100 --target 0 --nkin 200 --e0 10.5 --iFM 0 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_proton.lund"
        ;;
    copper)
        # target 15: Cu (Z=29, A=63), Accardi SVG supports 63Cu
        ARGS="--nevent 100 --target 15 --nkin 200 --e0 10.5 --iFM 2 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_copper.lund"
        ;;

    # ── iFM=5: R. Wiringa — supports 2H, 3H, 3He, 4He, 6Li, 7Li ─────
    deuterium)
        ARGS="--nevent 100 --target 1 --nkin 200 --e0 10.5 --iFM 5 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_deuterium.lund"
        ;;
    tritium)
        ARGS="--nevent 100 --target 2 --nkin 200 --e0 10.5 --iFM 5 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_tritium.lund"
        ;;
    he3)
        ARGS="--nevent 100 --target 3 --nkin 200 --e0 10.5 --iFM 5 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_he3.lund"
        ;;
    he4)
        ARGS="--nevent 100 --target 4 --nkin 200 --e0 10.5 --iFM 5 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_he4.lund"
        ;;
    li6)
        ARGS="--nevent 100 --target 5 --nkin 200 --e0 10.5 --iFM 5 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_li6.lund"
        ;;
    li7)
        ARGS="--nevent 100 --target 6 --nkin 200 --e0 10.5 --iFM 5 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_li7.lund"
        ;;

    # ── iFM=1: Fermi sphere + Bodek-Ritchie tail — any target with rFM>0
    carbon)
        ARGS="--nevent 100 --target 7 --nkin 200 --e0 10.5 --iFM 1 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_carbon.lund"
        ;;
    aluminum)
        ARGS="--nevent 100 --target 8 --nkin 200 --e0 10.5 --iFM 1 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_aluminum.lund"
        ;;
    iron)
        ARGS="--nevent 100 --target 9 --nkin 200 --e0 10.5 --iFM 1 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_iron.lund"
        ;;
    tin)
        ARGS="--nevent 100 --target 10 --nkin 200 --e0 10.5 --iFM 1 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_tin.lund"
        ;;
    lead)
        ARGS="--nevent 100 --target 11 --nkin 200 --e0 10.5 --iFM 1 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_lead.lund"
        ;;

    # ── iFM=2: Accardi SVG — supports Ne, Kr, Xe (and C, Al, Sn, Pb) ─
    neon)
        ARGS="--nevent 100 --target 12 --nkin 200 --e0 10.5 --iFM 2 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_neon.lund"
        ;;
    krypton)
        ARGS="--nevent 100 --target 13 --nkin 200 --e0 10.5 --iFM 2 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_krypton.lund"
        ;;
    xenon)
        ARGS="--nevent 100 --target 14 --nkin 200 --e0 10.5 --iFM 2 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_xenon.lund"
        ;;

    # ── iFM=3: Accardi CS — additional test for iron ─────────────────
    iron_cs)
        ARGS="--nevent 100 --target 9 --nkin 200 --e0 10.5 --iFM 3 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_iron_cs.lund"
        ;;

    # ── iFM=4: Hard sphere — additional test for carbon ──────────────
    carbon_hs)
        ARGS="--nevent 100 --target 7 --nkin 200 --e0 10.5 --iFM 4 --seed 42"
        GOLDEN="${GOLDEN_DIR}/ref_carbon_hs.lund"
        ;;

    *)
        echo "Unknown config: ${CONFIG}"
        echo "Available: proton, deuterium, tritium, he3, he4, li6, li7,"
        echo "           carbon, aluminum, iron, tin, lead, neon, krypton, xenon, copper,"
        echo "           iron_cs, carbon_hs"
        exit 1
        ;;
esac

OUTPUT="${TMP_DIR}/test_output.lund"

echo "=== Golden file test: ${CONFIG} ==="
echo "Running: ${SIM} ${ARGS} --output ${OUTPUT}"

# Run the simulation
${SIM} ${ARGS} --output "${OUTPUT}"

# Check golden file exists
if [ ! -f "${GOLDEN}" ]; then
    echo "Golden file not found: ${GOLDEN}"
    echo "Generate it with:"
    echo "  ${SIM} ${ARGS} --output ${GOLDEN}"
    exit 2
fi

# Try exact match first
if diff -q "${GOLDEN}" "${OUTPUT}" > /dev/null 2>&1; then
    echo "PASS: Exact match with golden file"
    exit 0
fi

# Fall back to tolerance-based comparison
echo "Exact match failed, trying tolerance-based comparison..."
COMPARE_SCRIPT="${SCRIPT_DIR}/compare_lund.py"
if [ -f "${COMPARE_SCRIPT}" ]; then
    python3 "${COMPARE_SCRIPT}" "${GOLDEN}" "${OUTPUT}"
    exit $?
else
    echo "FAIL: Output differs from golden file and compare_lund.py not found"
    echo "First difference:"
    diff "${GOLDEN}" "${OUTPUT}" | head -20
    exit 1
fi

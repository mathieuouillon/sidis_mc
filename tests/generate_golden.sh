#!/bin/bash
# Generate all golden reference files for regression tests.
# Usage: ./tests/generate_golden.sh [simulation_binary]
# Run from the project root directory.

set -euo pipefail

SIM="${1:-./build/simulation}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GOLDEN_DIR="${SCRIPT_DIR}/golden"

mkdir -p "${GOLDEN_DIR}"

echo "=== Generating golden reference files ==="
echo "Using: ${SIM}"
echo ""

# iFM=0: no Fermi motion
${SIM} --nevent 100 --target 0  --nkin 200 --e0 10.5 --iFM 0 --seed 42 --output "${GOLDEN_DIR}/ref_proton.lund"
echo "  [done] proton"

# iFM=5: R. Wiringa (2H, 3H, 3He, 4He, 6Li, 7Li)
${SIM} --nevent 100 --target 1  --nkin 200 --e0 10.5 --iFM 5 --seed 42 --output "${GOLDEN_DIR}/ref_deuterium.lund"
echo "  [done] deuterium"
${SIM} --nevent 100 --target 2  --nkin 200 --e0 10.5 --iFM 5 --seed 42 --output "${GOLDEN_DIR}/ref_tritium.lund"
echo "  [done] tritium"
${SIM} --nevent 100 --target 3  --nkin 200 --e0 10.5 --iFM 5 --seed 42 --output "${GOLDEN_DIR}/ref_he3.lund"
echo "  [done] he3"
${SIM} --nevent 100 --target 4  --nkin 200 --e0 10.5 --iFM 5 --seed 42 --output "${GOLDEN_DIR}/ref_he4.lund"
echo "  [done] he4"
${SIM} --nevent 100 --target 5  --nkin 200 --e0 10.5 --iFM 5 --seed 42 --output "${GOLDEN_DIR}/ref_li6.lund"
echo "  [done] li6"
${SIM} --nevent 100 --target 6  --nkin 200 --e0 10.5 --iFM 5 --seed 42 --output "${GOLDEN_DIR}/ref_li7.lund"
echo "  [done] li7"

# iFM=1: Fermi sphere + Bodek-Ritchie tail (C, Al, Fe, Sn, Pb)
${SIM} --nevent 100 --target 7  --nkin 200 --e0 10.5 --iFM 1 --seed 42 --output "${GOLDEN_DIR}/ref_carbon.lund"
echo "  [done] carbon"
${SIM} --nevent 100 --target 8  --nkin 200 --e0 10.5 --iFM 1 --seed 42 --output "${GOLDEN_DIR}/ref_aluminum.lund"
echo "  [done] aluminum"
${SIM} --nevent 100 --target 9  --nkin 200 --e0 10.5 --iFM 1 --seed 42 --output "${GOLDEN_DIR}/ref_iron.lund"
echo "  [done] iron"
${SIM} --nevent 100 --target 10 --nkin 200 --e0 10.5 --iFM 1 --seed 42 --output "${GOLDEN_DIR}/ref_tin.lund"
echo "  [done] tin"
${SIM} --nevent 100 --target 11 --nkin 200 --e0 10.5 --iFM 1 --seed 42 --output "${GOLDEN_DIR}/ref_lead.lund"
echo "  [done] lead"

# iFM=2: Accardi SVG (Ne, Kr, Xe, Cu)
${SIM} --nevent 100 --target 12 --nkin 200 --e0 10.5 --iFM 2 --seed 42 --output "${GOLDEN_DIR}/ref_neon.lund"
echo "  [done] neon"
${SIM} --nevent 100 --target 13 --nkin 200 --e0 10.5 --iFM 2 --seed 42 --output "${GOLDEN_DIR}/ref_krypton.lund"
echo "  [done] krypton"
${SIM} --nevent 100 --target 14 --nkin 200 --e0 10.5 --iFM 2 --seed 42 --output "${GOLDEN_DIR}/ref_xenon.lund"
echo "  [done] xenon"
${SIM} --nevent 100 --target 15 --nkin 200 --e0 10.5 --iFM 2 --seed 42 --output "${GOLDEN_DIR}/ref_copper.lund"
echo "  [done] copper"

# iFM=3: Accardi CS (iron)
${SIM} --nevent 100 --target 9  --nkin 200 --e0 10.5 --iFM 3 --seed 42 --output "${GOLDEN_DIR}/ref_iron_cs.lund"
echo "  [done] iron_cs"

# iFM=4: Hard sphere (carbon)
${SIM} --nevent 100 --target 7  --nkin 200 --e0 10.5 --iFM 4 --seed 42 --output "${GOLDEN_DIR}/ref_carbon_hs.lund"
echo "  [done] carbon_hs"

echo ""
echo "=== All 18 golden files generated in ${GOLDEN_DIR} ==="

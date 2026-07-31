# sidis_mc

**Semi-inclusive DIS (SIDIS) Monte-Carlo event generator for nuclear targets.** It generates electron-scattering events with **PYTHIA 6.4.28**, applies nuclear effects (Fermi motion, parton energy loss / "quenching", pT broadening, spectators), and writes **LUND**-format events — the input for the CLAS12 simulation chain (GEMC + coatjava reconstruction, driven by [`../background_merging`](../background_merging/README.md)).

The physics core is the original PYTHIA6 Fortran; everything around it (kinematics, Fermi motion, nuclear density, quenching, event processing, LUND writing) is modern C++20 in the `farm::` namespace.

## Build

Meson project (`simulation`), mixing **C / C++20 / Fortran**. Requires a Fortran compiler and **`libgfortran`** (linked explicitly for the PYTHIA6 object), built at `-O3`.

```bash
meson setup build
meson compile -C build          # → build/simulation
```

## Usage

```
simulation --nevent <N> --output <file> --target <0-15> --nkin <N> --e0 <GeV>
           [--iFM <0-5>] [--seed <N>] [--help]
```

| Flag | | Meaning |
|---|---|---|
| `--nevent`, `-n` | required | number of events to generate |
| `--output`, `-o` | required | LUND output file |
| `--target`, `-t` | required | target index, `0–15` (see below) |
| `--nkin` | required | events between kinematics/PYTHIA re-initialisations |
| `--e0` | required | beam (electron) energy in GeV |
| `--iFM` | default `5` | Fermi-motion model (see below) |
| `--seed`, `-s` | default time-based | RNG seed (seeds both the C++ RNG and PYTHIA) |

**Targets:** `0`=p `1`=2H `2`=3H `3`=He3 `4`=He4 `5`=Li6 `6`=Li7 `7`=C `8`=Al `9`=Fe `10`=Sn `11`=Pb `12`=Ne `13`=Kr `14`=Xe `15`=Cu

**Fermi motion (`--iFM`):** `0`=none `1`=Bodek–Ritchie `2`=Accardi-SVG `3`=Accardi-CS `4`=hard-sphere `5`=Wiringa

Example — 100k events on tin at 10.6 GeV:

```bash
./build/simulation -n 100000 -o sn_10gev.lund -t 10 --nkin 1000 --e0 10.6 -s 42
```

> The remaining physics knobs (`qhat` = 0.34 GeV²/fm, `ehat`, quenching weights SW vs Arleo, Woods–Saxon density, pT-broadening model, isospin/spectator switches) are **compile-time defaults** in the `SimConfig` struct in `src/main.cpp` — not CLI flags. Edit and rebuild to change them.

## Per-event pipeline

For each event (`Simulation::process_event`):

1. **PYTHIA** generates the DIS event (fragmentation deferred when quenching is on).
2. **Inverse Lorentz boost** back to the lab frame (when Fermi motion is enabled).
3. **Quenching** — for nuclear targets, an interaction point is sampled from the nuclear density and parton energy loss / pT broadening applied.
4. **Fragmentation** is then run in PYTHIA.
5. **Spectators** are added for light nuclei (targets 1–4, when enabled).
6. **DIS variables computed** and final-state particles extracted → written as a LUND event (header + one line per particle, with a sampled vertex `vz`).

## Layout

| Path | What |
|---|---|
| `src/main.cpp` | CLI, `SimConfig`, `Simulation` engine, `LundWriter` |
| `src/include/` | the C++ modules: `physics.h`, `pythia6.h`, `event_processor.h`, `fermi_motion.h`, `accardi_fm.h`, `nuclear_density.h`, `quenching.h` |
| `src/external/pythia6428.f` | PYTHIA 6.4.28 (the only Fortran) |
| `datafiles/fmacc/` | Accardi Fermi-motion tables (`fermimotion2.{CS,SVG}.tbl`) |
| `datafiles/fmrw/` | per-nucleus momentum distributions (`h2.momentum`, `h3.momentum`, …) |
| `datafiles/qweight/` | quenching-weight tables (`cont03.all`, `cont05.all`, …) |
| `scripts/plot_lund.py` | plot distributions from a LUND file |
| `tests/` | golden-file + statistical regression tests |

## Tests

```bash
meson test -C build
```

- **Golden-file tests** — 18 configurations (`proton`, `deuterium`, `tritium`, `he3`, `he4`, `li6`, `li7`, `carbon`, `aluminum`, `iron`, `tin`, `lead`, `neon`, `krypton`, `xenon`, `copper`, plus `iron_cs` and `carbon_hs` FM variants) run the generator with a fixed seed and diff the LUND output against stored references (`tests/test_golden.sh`, `tests/compare_lund.py`). Regenerate references with `tests/generate_golden.sh`.
- **Statistical distribution tests** — `tests/test_distributions.py` (needs `python3`).

⚠️ The LUND writer deliberately uses C `fprintf` with fixed format strings so output stays **byte-identical** to the golden files. Changing the formatting will break every golden test.

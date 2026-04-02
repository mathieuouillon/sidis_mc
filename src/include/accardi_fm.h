#ifndef FARM_ACCARDI_FM_H
#define FARM_ACCARDI_FM_H

// ============================================================================
// Accardi Fermi motion routines (iFM=2, iFM=3)
//
// C++ port of the rhofermi family of functions from fermimotion.f90
// Original author: Alberto Accardi (2005-2006)
//
// Models:
//   irho=1  SVG  (Sandel-Vary-Garpman, soft+hard)
//   irho=2  CS   (Ciofi-Simula parametrization)
//
// References:
//   [1] M. Sandel, J.P. Vary and S.I.A. Garpman, Phys.Rev.C20(1979)744
//   [2] C. Ciofi degli Atti and S. Simula, Phys.Rev.C53(1996)1689
//   [3] R.B. Wiringa et al., Phys.Rev.C62(2001)014001
// ============================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <span>
#include <stdexcept>
#include <sstream>
#include "physics.h"
#include <string>

namespace farm {

// ============================================================================
// Binary search utility (port of ilocatetab)
// ============================================================================

// Given an array xx[0..n-1] (monotonic) and a value x, returns j such that
// x is between xx[j] and xx[j+1].  j=-1 or j=n-1 indicates out of range.
// NOTE: Fortran used 1-based indexing; this version uses 0-based.
// The returned j is 0-based (j=0 means xx[0] matched, etc.)
// To match Fortran semantics exactly, we keep 1-based logic internally
// and convert at the boundaries.
inline void ilocatetab(std::span<const int> xx, int x, int& j) {
    int n = static_cast<int>(xx.size());
    // Fortran-style 1-based binary search
    int jl = 0;
    int ju = n + 1;
    while (ju - jl > 1) {
        int jm = (ju + jl) / 2;
        if ((xx[n - 1] >= xx[0]) == (x >= xx[jm - 1])) {
            jl = jm;
        } else {
            ju = jm;
        }
    }

    if (x == xx[0]) {
        j = 1;
    } else if (x == xx[n - 1]) {
        j = n - 1;
    } else {
        j = jl;
    }
    // j is now 1-based (0 means out of range on the low side)
}

// ============================================================================
// SVG soft part (reads tabulated parameters from file)
// ============================================================================

struct SVGSoftState {
    static constexpr int NMAX = 20;
    bool initialized = false;
    int nnuke = 0;
    int Zold = 0, Aold = 0;
    int n = 0; // current nucleus index (1-based)
    std::array<int, NMAX> ZZ{};
    std::array<int, NMAX> AA{};
    std::array<std::array<double, NMAX>, 6> a{};   // a[j][nuke], j=0..5, nuke=0..nnuke-1
    std::array<double, NMAX> alpha{};
    std::array<double, NMAX> beta{};
};

inline double rhosoftSVG(int iZ, int iA, double k, SVGSoftState& state) {
    constexpr double pi = farm::constants::pi;

    if (iA < 7) {
        throw std::runtime_error("ERROR (rhosoftSVG): called with A<7: " + std::to_string(iA)
                                 + " -- Try using Ciofi-Simula parametrization instead");
    }

    // Lazy initialization: read table file on first call
    if (!state.initialized) {
        state.initialized = true;

        std::ifstream fin("datafiles/fmacc/fermimotion2.SVG.tbl");
        if (!fin.is_open()) {
            throw std::runtime_error("ERROR: Cannot open datafiles/fmacc/fermimotion2.SVG.tbl");
        }

        // Skip header lines until column header containing "alpha"
        std::string line;
        while (std::getline(fin, line)) {
            if (line.find("alpha") != std::string::npos)
                break;
        }
        std::getline(fin, line); // skip blank line after column header

        // Read nucleus entries until sentinel (ZZ==0)
        int j = 0;
        while (j < SVGSoftState::NMAX) {
            std::string name;
            int zz, aa;
            double a1, a2, a3, a4, a5, a6, al, be;
            fin >> name >> zz >> aa >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> al >> be;
            if (zz == 0) break;
            state.ZZ[j] = zz;
            state.AA[j] = aa;
            state.a[0][j] = a1;
            state.a[1][j] = a2;
            state.a[2][j] = a3;
            state.a[3][j] = a4;
            state.a[4][j] = a5;
            state.a[5][j] = a6;
            state.alpha[j] = al;
            state.beta[j] = be;
            j++;
        }
        state.nnuke = j;
        fin.close();
    }

    // If new nucleus, locate nearest nucleus in table
    if (iZ != state.Zold || iA != state.Aold) {
        state.Zold = iZ;
        state.Aold = iA;
        int j;
        ilocatetab(std::span{state.ZZ.data(), static_cast<size_t>(state.nnuke)}, iZ, j);
        if (j == 0) {
            state.n = 1;
        } else if (state.ZZ[j - 1] == iZ && state.AA[j - 1] == iA) {
            // Exact match (j is 1-based, arrays are 0-based)
            state.n = j;
        } else {
            // Choose nucleus with closest A
            // j is 1-based; compare AA[j-1] and AA[j]
            if ((iA - state.AA[j - 1]) <= (state.AA[j] - iA)) {
                state.n = j;
            } else {
                state.n = j + 1;
            }
        }
    }

    // Compute the soft Fermi momentum distribution
    // n is 1-based; array index is n-1
    int idx = state.n - 1;
    double result = 0.0;
    for (int j = 0; j < 6; j++) {
        double nj = state.alpha[idx] * std::pow(state.beta[idx], j + 1);
        result += state.a[j][idx] * std::pow(nj / pi, 1.5) * std::exp(-nj * k * k);
    }

    return std::abs(result);
}

// ============================================================================
// SVG hard tail
// ============================================================================

struct SVGHardState {
    bool first_warned = false;
};

inline double rhohardSVG(int iZ, int iA, double k,
                          double& height, double& slope, double& norm,
                          SVGHardState& state) {
    constexpr double pi = farm::constants::pi;
    constexpr double fourpi = 4.0 * pi;

    if (!state.first_warned) {
        state.first_warned = true;
        if (iA < 8) {
            std::cerr << "**************************************************\n"
                  << "WARNING (rhohardSVG): called with 2<A<8: A=" << iA << "\n"
                  << " -- parametrization of hard tail's height\n"
                  << "    should be used with caution\n"
                  << "**************************************************\n";
        }
    }

    if (iZ == 2 && iA == 3) {
        slope = 0.234;
    } else {
        slope = 0.220;
    }

    if (iZ == 2 && iA == 4) {
        height = 0.0244 / fourpi;
    } else {
        double AA = static_cast<double>(iA);
        height = 0.00623 / fourpi
                 * (1.0 + 2.69 * std::pow(std::log(AA / 2.0), 0.695) / std::pow(AA, 0.133));
    }

    norm = 1.0 - height * std::pow(pi / slope, 1.5);

    return height * std::exp(-slope * k * k);
}

// ============================================================================
// SVG soft+hard front-end
// ============================================================================

inline double rhoSVGsh(int iZ, int iA, double k, int idist,
                        SVGSoftState& soft_state, SVGHardState& hard_state) {
    double norm, rs, rh, height, slope;

    if (idist == 1) {
        // soft+hard
        rs = rhosoftSVG(iZ, iA, k, soft_state);
        rh = rhohardSVG(iZ, iA, k, height, slope, norm, hard_state);
        return norm * rs + rh;
    } else if (idist == 2) {
        // soft (unnormalized)
        rs = rhosoftSVG(iZ, iA, k, soft_state);
        rh = rhohardSVG(iZ, iA, k, height, slope, norm, hard_state);
        return norm * rs;
    } else if (idist == 3) {
        // soft (normalized)
        return rhosoftSVG(iZ, iA, k, soft_state);
    } else if (idist == 4) {
        // hard (unnormalized)
        return rhohardSVG(iZ, iA, k, height, slope, norm, hard_state);
    }
    return 0.0;
}

// ============================================================================
// CS (Ciofi-Simula) parametrization
// ============================================================================

struct CSState {
    static constexpr int NMAX = 20;
    bool initialized = false;
    int nnuke = 0;
    int Zold = 0, Aold = 0;
    int n = 0; // current nucleus index (1-based)
    std::array<int, NMAX> ZZ{};
    std::array<int, NMAX> AA{};
    // n0(k) parameters
    std::array<double, NMAX> a0{}, b0{}, c0{};
    std::array<double, NMAX> d0{}, e0{}, f0{};
    // n1(k) parameters
    std::array<double, NMAX> a1{}, b1{}, b2{};
    std::array<double, NMAX> c1{}, d1{};
};

inline double rhoCSsh(int iZ, int iA, double k, int idist, CSState& state) {
    constexpr double pi = farm::constants::pi;
    constexpr double fourpi = 4.0 * pi;
    constexpr double pith = pi * 5.5683279968317;
    constexpr double normCS = 1.0 / (4.0 * pi);

    // Lazy initialization: read table file on first call
    if (!state.initialized) {
        state.initialized = true;

        std::ifstream fin("datafiles/fmacc/fermimotion2.CS.tbl");
        if (!fin.is_open()) {
            throw std::runtime_error("ERROR: Cannot open datafiles/fmacc/fermimotion2.CS.tbl");
        }

        // Skip header lines until we find the column header "Z   A"
        // then skip one more blank line before data starts
        std::string line;
        while (std::getline(fin, line)) {
            // Look for the column header line containing "Z" and "A"
            if (line.find("Z") != std::string::npos && line.find("a0") != std::string::npos)
                break;
        }
        std::getline(fin, line); // skip blank line after column header

        // Read n0(k) parameters until sentinel (ZZ==0)
        int j = 0;
        while (j < CSState::NMAX) {
            std::string name;
            int zz, aa;
            double va0, vb0, vc0, vd0, ve0, vf0;
            fin >> name >> zz >> aa >> va0 >> vb0 >> vc0 >> vd0 >> ve0 >> vf0;
            if (zz == 0) break;
            state.ZZ[j] = zz;
            state.AA[j] = aa;
            state.a0[j] = va0;
            state.b0[j] = vb0;
            state.c0[j] = vc0;
            state.d0[j] = vd0;
            state.e0[j] = ve0;
            state.f0[j] = vf0;
            j++;
        }
        state.nnuke = j;

        // Skip to n1 section: find column header with "a1"
        while (std::getline(fin, line)) {
            if (line.find("a1") != std::string::npos && line.find("b1") != std::string::npos)
                break;
        }
        std::getline(fin, line); // skip blank line after column header

        // Read n1(k) parameters
        j = 0;
        while (j < CSState::NMAX) {
            std::string name;
            int zz, aa;
            double va1, vb1, vb2, vc1, vd1;
            if (!(fin >> name >> zz >> aa >> va1 >> vb1 >> vb2 >> vc1 >> vd1)) break;
            if (zz == 0) break;
            // Verify ZZ/AA match (they should be the same nuclei)
            state.a1[j] = va1;
            state.b1[j] = vb1;
            state.b2[j] = vb2;
            state.c1[j] = vc1;
            state.d1[j] = vd1;
            j++;
        }
        // Don't overwrite nnuke — it was set from n0 section

        fin.close();
    }

    // If new nucleus, locate it in the tables
    if (iZ != state.Zold || iA != state.Aold) {
        state.Zold = iZ;
        state.Aold = iA;
        int j;
        ilocatetab(std::span{state.ZZ.data(), static_cast<size_t>(state.nnuke)}, iZ, j);
        // j is 1-based; find matching (Z,A) at j-1, j, or j+1 (Fortran-style)
        state.n = 0;
        for (int delta = 0; delta <= 1 && state.n == 0; delta++) {
            for (int sign = 0; sign <= 1 && state.n == 0; sign++) {
                int idx = j + (sign == 0 ? -delta : delta);
                if (idx >= 1 && idx <= state.nnuke) {
                    if (state.ZZ[idx - 1] == iZ && state.AA[idx - 1] == iA)
                        state.n = idx;
                }
            }
        }
        if (state.n == 0) {
            throw std::runtime_error("ERROR (rhoCSsh): Z,A outside parameter table: "
                                     + std::to_string(iZ) + " " + std::to_string(iA));
        }
    }

    // Compute the Fermi momentum distribution
    // n is 1-based; array index is n-1
    int idx = state.n - 1;
    double k2 = k * k;

    // n0 (as defined in [1])
    double n0;
    if (iA <= 4) {
        n0 = state.a0[idx] * std::exp(-state.b0[idx] * k2)
             / ((1.0 + state.c0[idx] * k2) * (1.0 + state.c0[idx] * k2))
             + state.d0[idx] * std::exp(-state.e0[idx] * k2)
             / ((1.0 + state.f0[idx] * k2) * (1.0 + state.f0[idx] * k2));
    } else {
        double k4 = k2 * k2;
        double k6 = k4 * k2;
        double k8 = k4 * k4;
        n0 = state.a0[idx] * std::exp(-state.b0[idx] * k2)
             * (1.0 + state.c0[idx] * k2 + state.d0[idx] * k4
                + state.e0[idx] * k6 + state.f0[idx] * k8);
    }

    // soft, hard, and n1
    double tmp1 = state.a1[idx] * std::exp(-state.b1[idx] * k2)
                  / ((1.0 + state.b2[idx] * k2) * (1.0 + state.b2[idx] * k2));
    double hard = state.c1[idx] * std::exp(-state.d1[idx] * k2);
    double soft = n0 + tmp1;
    double n1 = tmp1 + hard;

    // Final result
    if (idist == 1) {
        // soft+hard
        return normCS * (n0 + n1);
    } else if (idist == 2) {
        // soft (unnormalized)
        return normCS * soft;
    } else if (idist == 3) {
        // soft (normalized to 1)
        double norm_val = 1.0 - (state.c1[idx] / fourpi) * (pith / state.d1[idx]);
        return normCS * (soft / norm_val);
    } else if (idist == 4) {
        // hard (unnormalized)
        return normCS * hard;
    } else if (idist == 5) {
        // n0 [1]
        return normCS * n0;
    } else if (idist == 6) {
        // n1 [1]
        return normCS * n1;
    }
    return 0.0;
}

// CS convenience wrapper (soft+hard)
inline double rhoCS(int iZ, int iA, double k, CSState& state) {
    return rhoCSsh(iZ, iA, k, 1, state);
}

// ============================================================================
// Main dispatcher: rhofermi
// ============================================================================

struct AccardiFMState {
    SVGSoftState svg_soft;
    SVGHardState svg_hard;
    CSState cs;
};

// Returns the nucleon Fermi momentum distribution [GeV^-3]
//   iZ, iA  = atomic and mass numbers
//   k [GeV] = nucleon momentum
//   irho    = model: 1=SVG, 2=CS
inline double rhofermi(int iZ, int iA, double k, int irho,
                        AccardiFMState& state) {
    constexpr double FmGeV = farm::constants::hbar_c;
    constexpr double Fm3GeV3 = FmGeV * FmGeV * FmGeV;

    double kk = k / FmGeV; // convert GeV to fm^-1

    if (irho == 1) {
        return rhoSVGsh(iZ, iA, kk, 1, state.svg_soft, state.svg_hard) / Fm3GeV3;
    } else if (irho == 2) {
        return rhoCSsh(iZ, iA, kk, 1, state.cs) / Fm3GeV3;
    }
    return 0.0;
}

// ============================================================================
// GenFMtable: build cumulative FM table by integrating rhofermi
// ============================================================================

struct FMTableState {
    static constexpr int TABLE_SIZE = 1000;
    int FMnb = 0;
    float step_size = 0.0f;
    std::array<float, TABLE_SIZE> FM_table{};
};

// irho: 1=SVG, 2=CS
// iZ, iA: target nucleus
// iFM: Fermi motion model index (2 or 3, maps to irho = iFM-1)
// iTg: target index (used for CS substitutions)
// FMlimit: upper momentum cutoff [GeV]
inline void GenFMtable(int irho, int iZ, int iA, int iFM, int iTg,
                        float FMlimit, FMTableState& table,
                        AccardiFMState& fm_state) {
    constexpr double pi_val = farm::constants::pi;

    int itz = iZ;
    int ita = iA;

    // CS model substitutions for certain targets
    if (iTg == 12 && iFM == 3) {
        itz = 8;
        ita = 16;
    } else if (iFM == 3 && (iTg == 13 || iTg == 14)) {
        itz = 26;
        ita = 56;
    }

    table.FMnb = 1000;
    table.step_size = FMlimit / table.FMnb;
    double mom = 0.0;
    double ptot = 0.0;

    for (int i = 0; i < table.FMnb; i++) {
        double proba = pi_val * mom * mom
                       * rhofermi(itz, ita, mom, irho, fm_state)
                       * table.step_size;
        ptot += proba;
        table.FM_table[i] = static_cast<float>(ptot);
        mom += table.step_size;
    }

    // Normalize to 1
    float ptot_f = static_cast<float>(ptot);
    std::transform(table.FM_table.begin(), table.FM_table.begin() + table.FMnb,
                   table.FM_table.begin(),
                   [ptot_f](float v) { return v / ptot_f; });
}

} // namespace farm

#endif // FARM_ACCARDI_FM_H

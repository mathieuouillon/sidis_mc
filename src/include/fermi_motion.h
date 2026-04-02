#ifndef FARM_FERMI_MOTION_H
#define FARM_FERMI_MOTION_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <functional>
#include "physics.h"

namespace farm {

// ============================================================================
// Fermi motion tables and sampling
// ============================================================================

struct FermiMotionState {
    static constexpr int TABLE_SIZE = 1000;

    float step_size = 0.0f;
    float FMintact = 1.0f;
    std::array<float, TABLE_SIZE> FM_table{};
    std::array<float, TABLE_SIZE> FM_n{};
    std::array<float, TABLE_SIZE> FM_p{};
    std::array<float, TABLE_SIZE> FM_i{};
    int FMnb = 0;
};

// Read Wiringa momentum table from file
inline void gen_rw_table(FermiMotionState& fm, int iZ, int iA, float FMlimit) {
    constexpr float HBAR_C = constants::hbar_c_f;

    std::string filename;
    if (iA == 2 && iZ == 1)      filename = "datafiles/fmrw/h2.momentum";
    else if (iA == 3 && iZ == 1) filename = "datafiles/fmrw/h3.momentum";
    else if (iA == 3 && iZ == 2) filename = "datafiles/fmrw/he3.momentum";
    else if (iA == 4 && iZ == 2) filename = "datafiles/fmrw/he4.momentum";
    else if (iA == 6 && iZ == 3) filename = "datafiles/fmrw/lad.momentum";
    else if (iA == 7 && iZ == 3) filename = "datafiles/fmrw/lat.momentum";
    else {
        throw std::runtime_error("ERROR: Unsupported iZ=" + std::to_string(iZ) + " iA=" + std::to_string(iA) + " for RW table");
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("ERROR: Cannot open " + filename);
    }

    // Skip header until ** marker
    std::string line;
    while (std::getline(file, line)) {
        if (line.size() >= 3 && line[1] == '*' && line[2] == '*') break;
    }

    int i = 0;
    float sum_n = 0, sum_p = 0;
    float a = 0, b = 0, c = 0, d = 0, e = 0;

    if (iZ == 1 && iA == 2) {
        while (a < FMlimit / HBAR_C && std::getline(file, line)) {
            std::istringstream iss(line);
            if (!(iss >> a >> b >> c >> d)) break;
            fm.FM_n[i] = d;
            fm.FM_p[i] = d;
            sum_n += d;
            sum_p += d;
            i++;
        }
    } else if ((iZ == 1 && iA == 3) || (iZ == 2 && iA == 3) || (iZ == 3 && iA == 7)) {
        while (a < FMlimit / HBAR_C && std::getline(file, line)) {
            std::istringstream iss(line);
            if (!(iss >> a >> b >> c >> d >> e)) break;
            fm.FM_n[i] = d;
            fm.FM_p[i] = b;
            sum_n += d;
            sum_p += b;
            i++;
        }
    } else if (iZ == 3 && iA == 6) {
        while (a < FMlimit / HBAR_C && std::getline(file, line)) {
            std::istringstream iss(line);
            if (!(iss >> a >> b >> c)) break;
            fm.FM_n[i] = b;
            fm.FM_p[i] = b;
            sum_n += b;
            sum_p += b;
            i++;
        }
    } else if (iZ == 2 && iA == 4) {
        while (a < FMlimit / HBAR_C && std::getline(file, line)) {
            std::istringstream iss(line);
            if (!(iss >> a >> b >> c >> d >> e)) break;
            fm.FM_n[i] = b * (a + 0.05f) * (a + 0.05f);
            fm.FM_p[i] = b * (a + 0.05f) * (a + 0.05f);
            sum_n += fm.FM_n[i];
            sum_p += fm.FM_p[i];
            fm.FM_i[i] = d / b;
            i++;
        }
    }

    fm.step_size = FMlimit / static_cast<float>(i);

    // Normalize
    std::transform(fm.FM_n.begin(), fm.FM_n.begin() + i, fm.FM_n.begin(),
                   [sum_n](float v) { return v / sum_n; });
    std::transform(fm.FM_p.begin(), fm.FM_p.begin() + i, fm.FM_p.begin(),
                   [sum_p](float v) { return v / sum_p; });

    // Cumulative
    std::partial_sum(fm.FM_n.begin(), fm.FM_n.begin() + i, fm.FM_n.begin());
    std::partial_sum(fm.FM_p.begin(), fm.FM_p.begin() + i, fm.FM_p.begin());
}

// Sample Fermi momentum parameters
struct FMSample {
    float ThFM, PhiFM, Kf;
    float FMintact;
};

inline FMSample sample_fermi_motion(const FermiMotionState& fm, int iFM, float rFM,
                                     float FMlimit, int nucleon, int iTg,
                                     std::mt19937& rng) {
    std::uniform_real_distribution<float> uniform(0.0f, 1.0f);
    constexpr float PI = constants::pi_f;

    FMSample s;
    s.FMintact = 1.0f;
    s.ThFM = std::acos(2.0f * uniform(rng) - 1.0f);
    s.PhiFM = 2.0f * PI * uniform(rng);

    if (iFM == 4) {
        s.Kf = rFM * std::cbrt(uniform(rng));

    } else if (iFM == 1) {
        constexpr float a = 2.0f;
        float Thr = 1.0f - 6.0f * (rFM * a / PI) * (rFM * a / PI);
        float C = 4.0f / 3.0f * PI * rFM * rFM * rFM;
        float Ps = 5.0f;

        while (Ps > FMlimit || Ps < 0.0f) {
            s.Kf = uniform(rng);
            if (s.Kf <= Thr) {
                s.Kf = std::cbrt(3.0f * C * s.Kf / 4.0f / PI / Thr);
            } else {
                s.Kf = -rFM / ((s.Kf - Thr) * PI * C / 8.0f / (rFM * rFM * rFM * rFM * rFM) / (a * a) - 1.0f);
            }
            Ps = s.Kf;
        }

    } else if (iFM == 2 || iFM == 3) {
        float Rd = uniform(rng);
        int i = 0;
        while (fm.FM_table[i] < Rd) i++;
        s.Kf = static_cast<float>(i) * fm.step_size + uniform(rng) * fm.step_size;

    } else if (iFM == 5) {
        float Rd = uniform(rng);
        int i = 0;
        if (nucleon == 2112) {
            while (fm.FM_n[i] < Rd) i++;
        } else {
            while (fm.FM_p[i] < Rd) i++;
        }
        s.Kf = static_cast<float>(i) * fm.step_size + uniform(rng) * fm.step_size;
        if (iTg != 1) s.FMintact = fm.FM_i[i];

        if (s.Kf > FMlimit) {
            throw std::runtime_error("WARNING: Kf=" + std::to_string(s.Kf) + " > FMlimit=" + std::to_string(FMlimit));
        }
    } else {
        throw std::runtime_error("ERROR: Unsupported iFM=" + std::to_string(iFM));
    }

    return s;
}

// ============================================================================
// Lorentz boost + rotation for initial kinematics
// (replaces lorentz_fm in simulation_helpers.f90)
// ============================================================================

struct BoostParams {
    float BB1, B1x, B1y, B1z;
    float Thi, Phi;
};

inline BoostParams apply_initial_boost(Kinematics& k) {
    BoostParams bp;
    bp.BB1 = k.PPn / k.EEn;
    bp.B1x = k.Pnx / k.EEn;
    bp.B1y = k.Pny / k.EEn;
    bp.B1z = k.Pnz / k.EEn;

    lorentz_boost(k.EEe, k.PPe, k.Pex, k.Pey, k.Pez, bp.BB1, bp.B1x, bp.B1y, bp.B1z);
    lorentz_boost(k.EEn, k.PPn, k.Pnx, k.Pny, k.Pnz, bp.BB1, bp.B1x, bp.B1y, bp.B1z);

    // Rotate around Y
    bp.Thi = -std::atan2(k.Pex, k.Pez);
    float ct = std::cos(bp.Thi), st = std::sin(bp.Thi);
    // Electron
    float tmp = k.Pex * ct + k.Pez * st;
    k.Pez = -k.Pex * st + k.Pez * ct;
    k.Pex = tmp;
    // Nucleon
    tmp = k.Pnx * ct + k.Pnz * st;
    k.Pnz = -k.Pnx * st + k.Pnz * ct;
    k.Pnx = tmp;

    // Rotate around Z
    bp.Phi = -std::atan2(k.Pey, k.Pez);
    float cp = std::cos(bp.Phi), sp = std::sin(bp.Phi);
    // Electron
    tmp = k.Pey * cp + k.Pez * sp;
    k.Pez = -k.Pey * sp + k.Pez * cp;
    k.Pey = tmp;
    // Nucleon
    tmp = k.Pny * cp + k.Pnz * sp;
    k.Pnz = -k.Pny * sp + k.Pnz * cp;
    k.Pny = tmp;

    return bp;
}

// ============================================================================
// Full kinematics setup with retry loop
// (replaces farm_setup_kinematics)
// ============================================================================

struct KinematicsResult {
    BoostParams boost;
    int nucleon;
    float beam_energy;
    float nuc_theta, nuc_phi, nuc_momentum;
    float FMintact;
};

[[nodiscard]] inline KinematicsResult setup_kinematics(
    int ievent, int nevent, int iZ, int iA,
    int iFM, float rFM, float FMlimit, float E0,
    const FermiMotionState& fm, std::mt19937& rng)
{
    KinematicsResult result{};

    while (true) {
        result.nucleon = (ievent < nevent * iZ / iA) ? 2212 : 2112;

        float ThFM = 0, PhiFM = 0, Kf = 0;
        result.FMintact = 1.0f;
        if (rFM != 0.0f && iFM != 0) {
            auto s = sample_fermi_motion(fm, iFM, rFM, FMlimit, result.nucleon, 0, rng);
            ThFM = s.ThFM;
            PhiFM = s.PhiFM;
            Kf = s.Kf;
            result.FMintact = s.FMintact;
        }

        auto kin = init_kinematics(E0, Kf, ThFM, PhiFM);

        result.nuc_theta = ThFM;
        result.nuc_phi = PhiFM;
        result.nuc_momentum = Kf;

        if (iFM != 0) {
            result.boost = apply_initial_boost(kin);
        }

        result.beam_energy = kin.PPe;
        if (kin.PPe >= 4.0f) break;
    }

    return result;
}

} // namespace farm

#endif // FARM_FERMI_MOTION_H

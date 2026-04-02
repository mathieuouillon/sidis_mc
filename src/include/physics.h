#ifndef FARM_PHYSICS_H
#define FARM_PHYSICS_H

#include <cmath>
#include <numbers>
#include <string_view>

// ============================================================================
// Physical and mathematical constants
// ============================================================================

namespace farm::constants {
    inline constexpr double pi       = std::numbers::pi;
    inline constexpr float  pi_f     = std::numbers::pi_v<float>;
    inline constexpr double rad2deg  = 180.0 / std::numbers::pi;
    inline constexpr double deg2rad  = std::numbers::pi / 180.0;
    inline constexpr double hbar_c   = 0.1973269602;       // GeV*fm
    inline constexpr float  hbar_c_f = 0.1973269602f;
    inline constexpr float  electron_mass = 0.000511f;      // GeV
    inline constexpr float  nucleon_mass  = 0.938f;         // GeV
}

// ============================================================================
// Pure C++ physics functions
// ============================================================================

namespace farm {

// --- Lorentz transformation -------------------------------------------------
// Applies a Lorentz boost to a 4-vector (E, P, Px, Py, Pz) given
// boost velocity (BB, Bx, By, Bz) where BB = |beta|.
// The 4-vector is modified in place.
// NOTE: BB must be > 0 (division by BB^2).
inline void lorentz_boost(float& E, float& P, float& Px, float& Py, float& Pz,
                           float BB, float Bx, float By, float Bz) noexcept {
    float GG = 1.0f / std::sqrt(1.0f - BB * BB);

    float Ei = GG * (E - Bx * Px - By * Py - Bz * Pz);

    float Pix = -Bx * GG * E + Px + Px * (GG - 1) * Bx * Bx / (BB * BB)
                + Py * (GG - 1) * By * Bx / (BB * BB)
                + Pz * (GG - 1) * Bz * Bx / (BB * BB);

    float Piy = -By * GG * E + Px * (GG - 1) * Bx * By / (BB * BB)
                + Py + Py * (GG - 1) * By * By / (BB * BB)
                + Pz * (GG - 1) * Bz * By / (BB * BB);

    float Piz = -Bz * GG * E + Px * (GG - 1) * Bx * Bz / (BB * BB)
                + Pz + Py * (GG - 1) * By * Bz / (BB * BB)
                + Pz * (GG - 1) * Bz * Bz / (BB * BB);

    E  = Ei;
    Px = Pix;
    Py = Piy;
    Pz = Piz;
    P  = std::sqrt(Px * Px + Py * Py + Pz * Pz);
}

// --- Nuclear parameter lookup -----------------------------------------------
struct NuclearParams {
    int Z;         // atomic number
    int A;         // mass number
    float rFM;     // Fermi momentum parameter (GeV)
};

[[nodiscard]] constexpr NuclearParams get_nuclear_params(int target) {
    switch (target) {
        case 0:  return {1,  1,   0.0f};     // proton
        case 1:  return {1,  2,   0.07f};    // deuterium
        case 2:  return {1,  3,   0.1f};     // tritium
        case 3:  return {2,  3,   0.1f};     // 3He
        case 4:  return {2,  4,   0.120f};   // 4He
        case 5:  return {3,  6,   0.17f};    // 6Li
        case 6:  return {3,  7,   0.17f};    // 7Li
        case 7:  return {6,  12,  0.221f};   // Carbon
        case 8:  return {13, 27,  0.235f};   // Aluminum
        case 9:  return {26, 56,  0.260f};   // Iron
        case 10: return {50, 120, 0.260f};   // Tin
        case 11: return {82, 208, 0.265f};   // Lead
        case 12: return {10, 20,  0.228f};   // Neon
        case 13: return {36, 84,  0.260f};   // Krypton
        case 14: return {54, 132, 0.260f};   // Xenon
        case 15: return {29, 63,  0.260f};   // Copper
        default: return {1,  1,   0.0f};     // default (proton)
    }
}

// --- Vertex z position ------------------------------------------------------
// Generates random vertex z position for the given target type.
// Uses the provided random number (0-1) for position within target.
// For Carbon (dual target), uses random_choice (0-1) to select foil.
[[nodiscard]] constexpr float vertex_z(int target_type, float random_pos, float random_choice) noexcept {
    float target_pos, target_length;

    switch (target_type) {
        case 10: // Sn
            target_pos = -3.5f;
            target_length = 0.018f;
            break;
        case 15: // Cu
            target_pos = -8.5f;
            target_length = 0.009f;
            break;
        case 7:  // Carbon (dual foil)
            target_length = 0.2f;
            target_pos = (random_choice < 0.5f) ? -3.5f : -8.5f;
            break;
        case 1:  // LD2
            target_pos = -5.0f;
            target_length = 5.0f;
            break;
        default:
            return 0.0f;
    }
    return target_pos + target_length * (random_pos - 0.5f);
}

// --- Kinematics initialization ----------------------------------------------
struct Kinematics {
    float EEe, PPe, Pex, Pey, Pez;  // electron 4-momentum
    float EEn, PPn, Pnx, Pny, Pnz;  // nucleon 4-momentum
};

[[nodiscard]] inline Kinematics init_kinematics(float E0, float Kf, float ThFM, float PhiFM) noexcept {
    constexpr float electron_mass = constants::electron_mass;
    constexpr float nucleon_mass  = constants::nucleon_mass;

    Kinematics k;
    k.PPe = E0;
    k.EEe = std::sqrt(k.PPe * k.PPe + electron_mass * electron_mass);
    k.Pex = 0.0f;
    k.Pey = 0.0f;
    k.Pez = E0;

    k.PPn = Kf;
    k.EEn = std::sqrt(k.PPn * k.PPn + nucleon_mass * nucleon_mass);
    k.Pnx = std::sin(ThFM) * std::cos(PhiFM) * k.PPn;
    k.Pny = std::sin(ThFM) * std::sin(PhiFM) * k.PPn;
    k.Pnz = std::cos(ThFM) * k.PPn;

    return k;
}

} // namespace farm

#endif // FARM_PHYSICS_H

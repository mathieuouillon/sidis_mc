#ifndef FARM_EVENT_PROCESSOR_H
#define FARM_EVENT_PROCESSOR_H

#include <cmath>
#include <vector>
#include "pythia6.h"
#include "physics.h"

namespace farm {

// ============================================================================
// Particle data extracted from PYTHIA event record
// ============================================================================

struct Particle {
    int    id;          // PDG ID
    int    mother_id;   // PDG ID of mother
    int    charge;      // electric charge
    float  px, py, pz;  // 3-momentum
    float  p;           // |p|
    float  E;           // energy
    float  m;           // mass
    float  z;           // fragmentation variable
    float  theta;       // polar angle (degrees)
    float  phi;         // azimuthal angle (degrees)
    float  phih;        // hadron phi (degrees)
    float  tt;          // missing mass squared
    float  Pts;         // transverse momentum squared
    float  Xf;          // Feynman x
    float  alpha_s;     // Nachtmann variable
};

struct EventData {
    float Nu = 0, Q2 = 0, xBj = 0, W = 0, y = 0;
    int   trk_gs = 0;   // virtual photon index
    std::vector<Particle> particles;
};

// ============================================================================
// Event processing: compute DIS variables and extract particles
// (replaces ComputV in book.f90)
// ============================================================================

inline EventData compute_event(Pythia6& py) {
    EventData evt;

    double eip  = py.get_p(1, 4);   // initial electron energy
    double nip  = 0.0;               // nucleon initial momentum (fixed target)
    double nie  = py.get_p(2, 5);    // target nucleon mass

    int N = py.n_particles();

    // First pass: find virtual photon and compute DIS variables
    for (int ip = 1; ip <= N; ip++) {
        if (py.pdg_id(ip) == 22 && py.mother(ip) == 1) {
            double p4 = py.get_p(ip, 4);
            double p3 = py.get_p(ip, 3);
            double p2 = py.get_p(ip, 2);
            double p1 = py.get_p(ip, 1);
            double p5 = py.get_p(ip, 5);

            evt.Nu  = static_cast<float>((p4 * nie + p3 * nip) / py.get_p(2, 5));
            evt.Q2  = static_cast<float>(p5 * p5);
            if (evt.Nu != 0.0f) {
                evt.xBj = evt.Q2 / (2.0f * static_cast<float>(py.get_p(2, 5)) * evt.Nu);
                evt.y   = evt.Nu * static_cast<float>(py.get_p(2, 5)) / static_cast<float>(eip * (nie + nip));
            }
            evt.W = static_cast<float>(std::sqrt(
                (nie + p4) * (nie + p4) - (-nip + p3) * (-nip + p3) - p2 * p2 - p1 * p1));
            evt.trk_gs = ip;
        }
    }

    if (evt.trk_gs == 0) return evt;  // no virtual photon found

    double phi_ele = std::atan2(py.get_p(evt.trk_gs, 2), py.get_p(evt.trk_gs, 1)) * 57.2958 + 210.0;

    // Second pass: extract final-state particles
    for (int ip = 1; ip <= N; ip++) {
        int kid = py.pdg_id(ip);
        int kst = py.status(ip);
        int abskid = std::abs(kid);

        bool select = (abskid == 211)
            || (kid == 11 && kst == 1)
            || (abskid == 321) || (kid == 310)
            || (abskid == 2212 && kst == 1)
            || (abskid == 2112 && kst == 1)
            || (abskid > 10000 && kst == 1)
            || (kid == 22) || (kid == 111);

        if (!select) continue;

        Particle part{};
        // Charge
        if (kid == 11 || kid == -211 || kid == -321 || kid == -2212)
            part.charge = -1;
        else if (kid == 211 || kid == 321 || kid == 1000010030 || kid == 1000010020 || kid == 2212)
            part.charge = 1;
        else if (kid == 1000020030)
            part.charge = 2;
        else
            part.charge = 0;

        double p1 = py.get_p(ip, 1);
        double p2 = py.get_p(ip, 2);
        double p3 = py.get_p(ip, 3);
        double p4 = py.get_p(ip, 4);
        double p5 = py.get_p(ip, 5);

        part.id        = kid;
        part.mother_id = py.pdg_id(py.mother(ip));
        part.px = static_cast<float>(p1);
        part.py = static_cast<float>(p2);
        part.pz = static_cast<float>(p3);
        part.p  = static_cast<float>(std::sqrt(p4 * p4 - p5 * p5));
        part.E  = static_cast<float>(p4);
        part.m  = static_cast<float>(p5);

        part.z = static_cast<float>((p4 * nie + p3 * nip) / (evt.Nu * py.get_p(2, 5)));

        double pmag = std::sqrt(p1 * p1 + p2 * p2 + p3 * p3);
        part.theta = (pmag > 0) ? static_cast<float>(57.2957795 * std::acos(p3 / pmag)) : 0.0f;
        part.phi = static_cast<float>(std::atan2(p2, p1) * 57.2958 + 30.0);
        if (part.phi < 0) part.phi += 360.0f;

        // Phih calculation
        double A1 = std::sin(phi_ele / 57.2958), A2 = -std::cos(phi_ele / 57.2958), A3 = 0;
        double AA = A1 * A1 + A2 * A2;
        double gp1 = py.get_p(evt.trk_gs, 1), gp2 = py.get_p(evt.trk_gs, 2), gp3 = py.get_p(evt.trk_gs, 3);
        double B1 = gp2 * p3 - gp3 * p2;
        double B2 = gp3 * p1 - gp1 * p3;
        double B3 = gp1 * p2 - gp2 * p1;
        double BB = B1 * B1 + B2 * B2 + B3 * B3;
        part.phih = (AA * BB > 0) ? static_cast<float>(std::acos((A1 * B1 + A2 * B2 + A3 * B3) / std::sqrt(AA * BB)) * 57.2958) : 0.0f;

        // tt (missing mass squared)
        double gp4 = py.get_p(evt.trk_gs, 4);
        part.tt = static_cast<float>((evt.Nu - p4) * (evt.Nu - p4)
            - (gp1 - p1) * (gp1 - p1) - (gp2 - p2) * (gp2 - p2) - (gp3 - p3) * (gp3 - p3));

        // Pts
        double zNu = part.z * evt.Nu;
        double dot4 = -gp4 * p4 + gp1 * p1 + gp2 * p2 + gp3 * p3;
        part.Pts = static_cast<float>(zNu * zNu - p5 * p5
            - (dot4 + zNu * evt.Nu) * (dot4 + zNu * evt.Nu) / (evt.Nu * evt.Nu + evt.Q2));

        // Xf
        double NuQ2 = evt.Nu * evt.Nu + evt.Q2;
        double Whalf = evt.W / 2.0;
        double denom = std::sqrt(Whalf * Whalf - p5 * p5) * evt.W * std::sqrt(NuQ2);
        if (std::abs(denom) > 1e-10) {
            part.Xf = static_cast<float>(
                (part.z * py.get_p(2, 5) * evt.Nu * evt.Nu
                 - part.z * evt.Q2 * evt.Nu
                 - (py.get_p(2, 5) + evt.Nu) * (-dot4)) / denom);
        }

        // Alpha_s (Nachtmann)
        double denom2 = std::sqrt(Whalf * Whalf + evt.Q2);
        if (std::abs(denom2) > 1e-10 && std::abs(p5) > 1e-10) {
            part.alpha_s = static_cast<float>(
                (p4 - (zNu * evt.Nu + dot4) / denom2) / p5);
        }

        evt.particles.push_back(part);
    }

    return evt;
}

// ============================================================================
// Rotate all PYTHIA particles around Y axis (replaces FinalRotY)
// ============================================================================

inline void rotate_all_y(Pythia6& py, float theta) {
    float ct = std::cos(theta), st = std::sin(theta);
    for (int ip = 1; ip <= py.n_particles(); ip++) {
        double p1 = py.get_p(ip, 1);
        double p3 = py.get_p(ip, 3);
        py.set_p(ip, 1, p1 * ct + p3 * st);
        py.set_p(ip, 3, -p1 * st + p3 * ct);
    }
}

// ============================================================================
// Rotate all PYTHIA particles around Z axis (replaces FinalRotZ)
// ============================================================================

inline void rotate_all_z(Pythia6& py, float phi) {
    float cp = std::cos(phi), sp = std::sin(phi);
    for (int ip = 1; ip <= py.n_particles(); ip++) {
        double p2 = py.get_p(ip, 2);
        double p3 = py.get_p(ip, 3);
        py.set_p(ip, 2, p2 * cp + p3 * sp);
        py.set_p(ip, 3, -p2 * sp + p3 * cp);
    }
}

// ============================================================================
// Inverse Lorentz boost all PYTHIA particles (replaces lorentz_fm_back)
// ============================================================================

inline void boost_all_back(Pythia6& py, float BB1, float B1x, float B1y, float B1z,
                            float Thi, float Phi) {
    rotate_all_z(py, -Phi);
    rotate_all_y(py, -Thi);

    for (int ip = 1; ip <= py.n_particles(); ip++) {
        float mom1 = static_cast<float>(py.get_p(ip, 1));
        float mom2 = static_cast<float>(py.get_p(ip, 2));
        float mom3 = static_cast<float>(py.get_p(ip, 3));
        float mom4 = static_cast<float>(py.get_p(ip, 4));
        float ppp;
        lorentz_boost(mom4, ppp, mom1, mom2, mom3, -BB1, -B1x, -B1y, -B1z);
        py.set_p(ip, 1, mom1);
        py.set_p(ip, 2, mom2);
        py.set_p(ip, 3, mom3);
        py.set_p(ip, 4, mom4);
    }
}

// ============================================================================
// Add nuclear spectator (replaces create_spec)
// ============================================================================

inline void add_spectator(Pythia6& py, int iTg, int nucleon,
                           float nuc_the, float nuc_phi, float nuc_mom,
                           float FMintact, float rand_val) {
    if (rand_val > FMintact) return;

    int specId = 0;
    switch (iTg) {
        case 1: specId = (nucleon == 2112) ? 2212 : 2112; break;
        case 2: specId = (nucleon == 2112) ? 1000010020 : 1000000020; break;
        case 3: specId = (nucleon == 2112) ? 1000020020 : 1000010020; break;
        case 4: specId = (nucleon == 2112) ? 1000020030 : 1000010030; break;
        default: return;
    }

    double spx = -std::sin(nuc_the) * std::cos(nuc_phi) * nuc_mom;
    double spy = -std::sin(nuc_the) * std::sin(nuc_phi) * nuc_mom;
    double spz = -std::cos(nuc_the) * nuc_mom;

    double mass = 0.0;
    switch (specId) {
        case 2212:       mass = 0.938272; break;
        case 2112:       mass = 0.939566; break;
        case 1000000020: mass = 1.87913;  break;
        case 1000010020: mass = 1.876124; break;
        case 1000020020: mass = 1.87654;  break;
        case 1000020030: mass = 2.809356; break;
        case 1000010030: mass = 2.809356; break;
    }

    double e = std::sqrt(spx * spx + spy * spy + spz * spz + mass * mass);
    py.add_particle(1, specId, 2, spx, spy, spz, e, mass);
}

} // namespace farm

#endif // FARM_EVENT_PROCESSOR_H

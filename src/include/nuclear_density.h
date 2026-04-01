#ifndef FARM_NUCLEAR_DENSITY_H
#define FARM_NUCLEAR_DENSITY_H

// C++ port of src/physics/nucdens.f90 and src/physics/density.f90
// Original Fortran by Alberto Accardi and Daniel Gruenewald

#include <cmath>
#include <functional>
#include <random>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <cstdio>

namespace farm {

// =====================================================================
// Woods-Saxon module state (replaces Fortran module nucdens_data_module)
// =====================================================================
struct WoodsSaxonParams {
    double RR = 0.0;   // radius parameter
    double a0 = 0.0;   // diffuseness
    double c0 = 0.0;   // deformation
};

// =====================================================================
// Adaptive Gaussian Quadrature (DGAUSS11)
// =====================================================================
template <typename F>
inline double dgauss11(F func, double A, double B, double EPS) {
    static constexpr double W[12] = {
        0.1012285362903762590e0,
        0.2223810344533744710e0,
        0.3137066458778872870e0,
        0.3626837833783619830e0,
        0.2715245941175409490e-1,
        0.6225352393864789290e-1,
        0.9515851168249278480e-1,
        0.1246289712555338720e0,
        0.1495959888165767320e0,
        0.1691565193950025380e0,
        0.1826034150449235890e0,
        0.1894506104550684960e0
    };
    static constexpr double X[12] = {
        0.9602898564975362320e0,
        0.7966664774136267400e0,
        0.5255324099163289860e0,
        0.1834346424956498050e0,
        0.9894009349916499330e0,
        0.9445750230732325760e0,
        0.8656312023878317440e0,
        0.7554044083550030340e0,
        0.6178762444026437480e0,
        0.4580167776572273860e0,
        0.2816035507792589130e0,
        0.9501250983763744020e-1
    };

    double result = 0.0;
    if (B == A) return result;
    double CONST_ = 0.005 / (B - A);
    double BB = A;

    // Outer loop: advance through intervals
    for (;;) {
        double AA = BB;
        BB = B;
        // Inner loop: subdivide until convergence
        for (;;) {
            double C1 = 0.5 * (BB + AA);
            double C2 = 0.5 * (BB - AA);
            double S8 = 0.0;
            for (int I = 0; I < 4; ++I) {
                double U = C2 * X[I];
                S8 += W[I] * (func(C1 + U) + func(C1 - U));
            }
            S8 *= C2;
            double S16 = 0.0;
            for (int I = 4; I < 12; ++I) {
                double U = C2 * X[I];
                S16 += W[I] * (func(C1 + U) + func(C1 - U));
            }
            S16 *= C2;
            if (std::abs(S16 - S8) <= EPS * (1.0 + std::abs(S16))) {
                // Converged
                result += S16;
                if (BB != B) {
                    break; // cycle outer
                } else {
                    return result;
                }
            }
            BB = C1;
            if (1.0 + std::abs(CONST_ * C2) == 1.0) {
                // Accuracy not achievable
                std::fprintf(stderr,
                    "    FUNCTION DGAUSS11 ... TOO HIGH ACCURACY REQUIRED\n");
                return 0.0;
            }
        }
    }
}

// =====================================================================
// Binary search in ordered table (locatetable)
// =====================================================================
// Given an array xx[0..n-1] (monotonic) and a value x, returns j such
// that x is between xx[j] and xx[j+1].  Uses 0-based indexing.
inline int locatetable(const double* xx, int n, double x) {
    int jl = -1;
    int ju = n;
    while (ju - jl > 1) {
        int jm = (ju + jl) / 2;
        if ((xx[n - 1] >= xx[0]) == (x >= xx[jm])) {
            jl = jm;
        } else {
            ju = jm;
        }
    }
    if (x == xx[0]) {
        return 0;
    } else if (x == xx[n - 1]) {
        return n - 2;
    } else {
        return jl;
    }
}

// =====================================================================
// Reid soft-core deuteron wave function [fm^-3]
// REF: R.V.Reid, Ann.Phys.(NY)50(68)411-448
// =====================================================================
inline double reidsc(double r) {
    static constexpr double fourpi = 2.0 * 6.2831530718;
    static constexpr double mu = 0.7;
    static constexpr double mu2 = mu * mu;
    static constexpr double AS2 = 0.7749985;
    static constexpr double ADS2 = 6.7081e-4;
    static constexpr double alpha = 0.33088;
    static constexpr double twoalpha = 0.66176;

    static constexpr int N = 33;

    static constexpr double xx[N] = {
        1.000e-2, 4.125e-2, 7.250e-2, 1.350e-1, 1.975e-1,
        2.600e-1, 3.225e-1, 3.850e-1, 4.475e-1, 5.100e-1,
        5.725e-1, 6.350e-1, 6.975e-1, 7.600e-1, 8.850e-1,
        1.010e0 , 1.135e0 , 1.260e0 , 1.385e0 , 1.510e0 ,
        1.760e0 , 2.010e0 , 2.510e0 , 3.010e0 , 3.510e0 ,
        4.010e0 , 4.510e0 , 5.010e0 , 5.510e0 , 6.010e0 ,
        7.010e0 , 8.010e0 , 9.010e0
    };
    static constexpr double u[N] = {
        0.0000e0, 3.3373e-5, 2.3901e-4, 2.7621e-3, 1.2737e-2,
        3.6062e-2, 7.5359e-2, 1.2847e-1, 1.8993e-1, 2.5349e-1,
        3.1390e-1, 3.6770e-1, 4.1317e-1, 4.4992e-1, 4.9953e-1,
        5.2406e-1, 5.3166e-1, 5.2846e-1, 5.1926e-1, 5.0621e-1,
        4.7505e-1, 4.4200e-1, 3.7864e-1, 3.2249e-1, 2.7399e-1,
        2.3251e-1, 1.9719e-1, 1.6718e-1, 1.4172e-1, 1.2012e-1,
        8.6290e-2, 6.1983e-2, 4.4523e-2
    };
    static constexpr double w[N] = {
        0.0000e0 , 1.0850e-5, 8.4073e-5, 1.0369e-3, 4.9642e-3,
        1.4446e-2, 3.0795e-2, 5.3157e-2, 7.8995e-2, 1.0525e-1,
        1.2933e-1, 1.4958e-1, 1.6529e-1, 1.7645e-1, 1.8710e-1,
        1.8654e-1, 1.7946e-1, 1.6910e-1, 1.5742e-1, 1.4553e-1,
        1.2314e-1, 1.0373e-1, 7.3859e-2, 5.3293e-2, 3.9077e-2,
        2.9115e-2, 2.2016e-2, 1.6871e-2, 1.3079e-2, 1.0243e-2,
        6.4412e-3, 4.1575e-3, 2.7363e-3
    };
    static constexpr double du[N] = {
        2.9751e-4, 2.6127e-3, 1.2335e-2, 8.2951e-2, 2.5326e-1,
        5.0052e-1, 7.5072e-1, 9.3349e-1, 1.0162e0 , 1.0034e0 ,
        9.2042e-1, 7.9674e-1, 6.5744e-1, 5.1985e-1, 2.8494e-1,
        1.1865e-1, 1.1397e-2,-5.4090e-2,-9.2523e-2,-1.1418e-1,
       -1.3097e-1,-1.3193e-1,-1.2004e-1,-1.0453e-1,-8.9709e-2,
       -7.6510e-2,-6.5054e-2,-5.5229e-2,-4.6850e-2,-3.9726e-2,
       -2.8547e-2,-2.0508e-2,-1.4731e-2
    };
    static constexpr double dw[N] = {
        7.4202e-5, 8.9321e-4, 4.4755e-3, 3.1892e-2, 1.0121e-1,
        2.0596e-1, 3.1477e-1, 3.9371e-1, 4.2466e-1, 4.0842e-1,
        3.5772e-1, 2.8848e-1, 2.1428e-1, 1.4427e-1, 3.3294e-2,
       -3.5899e-2,-7.3151e-2,-9.0091e-2,-9.5299e-2,-9.4158e-2,
       -8.3973e-2,-7.1282e-2,-4.9327e-2,-3.3940e-2,-2.3612e-2,
       -1.6691e-2,-1.2002e-2,-8.7768e-3,-6.5201e-3,-4.9136e-3,
       -2.8956e-3,-1.7758e-3,-1.1227e-3
    };

    double xval = mu * r;
    double result;

    if (xval < xx[1]) {
        result = (u[1] * u[1] + w[1] * w[1]);
    } else if (xval > xx[N - 1]) {
        double alphax = alpha * xval;
        double term = 1.0 + ADS2 * (1.0 + (3.0 / alphax) + (3.0 / (alphax * alphax)));
        result = AS2 * std::exp(-twoalpha * xval) * (term * term);
    } else {
        int j = locatetable(xx, N, xval);
        double h = xx[j + 1] - xx[j];
        double p = (xval - xx[j]) / h;
        double A1 = u[j];
        double A2 = h * du[j];
        double A3 = 3.0 * (u[j + 1] - u[j]) - (2.0 * du[j] + du[j + 1]) * h;
        double A4 = 2.0 * (u[j] - u[j + 1]) + (du[j] + du[j + 1]) * h;
        double u2_val = A1 + A2 * p + A3 * p * p + A4 * p * p * p;
        u2_val *= u2_val;

        A1 = w[j];
        A2 = h * dw[j];
        A3 = 3.0 * (w[j + 1] - w[j]) - (2.0 * dw[j] + dw[j + 1]) * h;
        A4 = 2.0 * (w[j] - w[j + 1]) + (dw[j] + dw[j + 1]) * h;
        double w2_val = A1 + A2 * p + A3 * p * p + A4 * p * p * p;
        w2_val *= w2_val;

        result = u2_val + w2_val;
    }

    return result;
}

// =====================================================================
// Deuteron density (based on Reid soft-core potential)
// =====================================================================
inline double DeuteronDensity(double r) {
    static constexpr double twopi = 6.2831530718;
    return reidsc(2.0 * r) / (twopi * r * r);
}

// r^2 * DeuteronDensity helper for integration
inline double r2DD(double r) {
    return 2.0 * r * r * reidsc(2.0 * r);
}

// =====================================================================
// Woods-Saxon 3D distribution (unnormalized, with 4*pi*r^2 Jacobian)
// Requires external WoodsSaxonParams for RR, a0, c0 state.
// =====================================================================
inline double WoodsSaxon3d(double r, const WoodsSaxonParams& ws) {
    static constexpr double pi = 3.141592653;
    double val = 4.0 * pi * r * r
        * (1.0 + ws.c0 * (r * r / (ws.RR * ws.RR)))
        / (1.0 + std::exp((r - ws.RR) / ws.a0));
    if (val < 0.0) val = 0.0;
    return val;
}

// r^2 * WoodsSaxon3d helper for integration
inline double r2WS3d(double r, const WoodsSaxonParams& ws) {
    return r * r * WoodsSaxon3d(r, ws);
}

// =====================================================================
// Woods-Saxon parameter tables (experimental values, Z=2..92)
// DATA statements from Fortran preserved exactly.
// =====================================================================
struct WSParamEntry {
    double R;   // radius parameter (fm)
    double a;   // diffuseness (fm)
    double c;   // deformation parameter
};

// Returns experimental WS parameters for a given Z (2..92).
// If no experimental data available (all zeros), returns {0,0,0}.
inline WSParamEntry get_ws_params(int Z) {
    // Table indexed by Z.  Only entries with known experimental values
    // are non-zero.  This reproduces the Fortran DATA statements exactly.
    static constexpr WSParamEntry table[93] = {
        /* Z= 0 */ {0.0, 0.0, 0.0},
        /* Z= 1 */ {0.0, 0.0, 0.0},
        /* Z= 2 */ {1.01,  0.327, 0.445},
        /* Z= 3 */ {0.0, 0.0, 0.0},
        /* Z= 4 */ {0.0, 0.0, 0.0},
        /* Z= 5 */ {0.0, 0.0, 0.0},
        /* Z= 6 */ {0.0, 0.0, 0.0},
        /* Z= 7 */ {2.570, 0.505, -0.180},
        /* Z= 8 */ {2.608, 0.513, -0.051},
        /* Z= 9 */ {2.58,  0.567, 0.0},
        /* Z=10 */ {2.791, 0.698, -0.168},
        /* Z=11 */ {0.0, 0.0, 0.0},
        /* Z=12 */ {3.192, 0.604, -0.249},
        /* Z=13 */ {3.07,  0.519, 0.0},
        /* Z=14 */ {3.340, 0.580, -0.233},
        /* Z=15 */ {3.369, 0.582, -0.173},
        /* Z=16 */ {3.20,  0.59,  0.0},
        /* Z=17 */ {3.476, 0.599, -0.10},
        /* Z=18 */ {3.73,  0.63,  -0.19},
        /* Z=19 */ {3.743, 0.585, -0.201},
        /* Z=20 */ {3.766, 0.586, -0.161},
        /* Z=21 */ {0.0, 0.0, 0.0},
        /* Z=22 */ {3.843, 0.588, 0.0},
        /* Z=23 */ {3.94,  0.505, 0.0},
        /* Z=24 */ {3.98,  0.542, 0.0},
        /* Z=25 */ {3.89,  0.567, 0.0},
        /* Z=26 */ {4.111, 0.558, 0.0},
        /* Z=27 */ {4.158, 0.575, 0.0},
        /* Z=28 */ {4.309, 0.517, -0.131},
        /* Z=29 */ {4.218, 0.596, 0.0},
        /* Z=30 */ {4.285, 0.584, 0.0},
        /* Z=31 */ {0.0, 0.0, 0.0},
        /* Z=32 */ {4.45,  0.573, 0.0},
        /* Z=33 */ {0.0, 0.0, 0.0},
        /* Z=34 */ {0.0, 0.0, 0.0},
        /* Z=35 */ {0.0, 0.0, 0.0},
        /* Z=36 */ {0.0, 0.0, 0.0},
        /* Z=37 */ {0.0, 0.0, 0.0},
        /* Z=38 */ {4.83,  0.496, 0.0},
        /* Z=39 */ {4.76,  0.571, 0.0},
        /* Z=40 */ {0.0, 0.0, 0.0},
        /* Z=41 */ {4.87,  0.573, 0.0},
        /* Z=42 */ {0.0, 0.0, 0.0},
        /* Z=43 */ {0.0, 0.0, 0.0},
        /* Z=44 */ {0.0, 0.0, 0.0},
        /* Z=45 */ {0.0, 0.0, 0.0},
        /* Z=46 */ {0.0, 0.0, 0.0},
        /* Z=47 */ {0.0, 0.0, 0.0},
        /* Z=48 */ {5.38,  0.532, 0.0},
        /* Z=49 */ {5.357, 0.563, 0.0},
        /* Z=50 */ {5.442, 0.543, 0.0},
        /* Z=51 */ {5.32,  0.57,  0.0},
        /* Z=52 */ {0.0, 0.0, 0.0},
        /* Z=53 */ {0.0, 0.0, 0.0},
        /* Z=54 */ {0.0, 0.0, 0.0},
        /* Z=55 */ {0.0, 0.0, 0.0},
        /* Z=56 */ {0.0, 0.0, 0.0},
        /* Z=57 */ {5.71,  0.535, 0.0},
        /* Z=58 */ {0.0, 0.0, 0.0},
        /* Z=59 */ {0.0, 0.0, 0.0},
        /* Z=60 */ {5.626, 0.618, 0.0},
        /* Z=61 */ {0.0, 0.0, 0.0},
        /* Z=62 */ {5.804, 0.581, 0.0},
        /* Z=63 */ {0.0, 0.0, 0.0},
        /* Z=64 */ {5.930, 0.576, 0.0},
        /* Z=65 */ {0.0, 0.0, 0.0},
        /* Z=66 */ {0.0, 0.0, 0.0},
        /* Z=67 */ {6.18,  0.57,  0.0},
        /* Z=68 */ {5.98,  0.446, 0.19},
        /* Z=69 */ {0.0, 0.0, 0.0},
        /* Z=70 */ {6.127, 0.363, 0.0},
        /* Z=71 */ {0.0, 0.0, 0.0},
        /* Z=72 */ {0.0, 0.0, 0.0},
        /* Z=73 */ {6.38,  0.64,  0.0},
        /* Z=74 */ {6.51,  0.535, 0.0},
        /* Z=75 */ {0.0, 0.0, 0.0},
        /* Z=76 */ {0.0, 0.0, 0.0},
        /* Z=77 */ {0.0, 0.0, 0.0},
        /* Z=78 */ {0.0, 0.0, 0.0},
        /* Z=79 */ {6.38,  0.535, 0.0},
        /* Z=80 */ {0.0, 0.0, 0.0},
        /* Z=81 */ {0.0, 0.0, 0.0},
        /* Z=82 */ {6.62,  0.546, 0.0},
        /* Z=83 */ {6.75,  0.468, 0.0},
        /* Z=84 */ {0.0, 0.0, 0.0},
        /* Z=85 */ {0.0, 0.0, 0.0},
        /* Z=86 */ {0.0, 0.0, 0.0},
        /* Z=87 */ {0.0, 0.0, 0.0},
        /* Z=88 */ {0.0, 0.0, 0.0},
        /* Z=89 */ {0.0, 0.0, 0.0},
        /* Z=90 */ {6.792, 0.571, 0.0},
        /* Z=91 */ {0.0, 0.0, 0.0},
        /* Z=92 */ {6.874, 0.556, 0.0}
    };
    if (Z < 0 || Z > 92) return {0.0, 0.0, 0.0};
    return table[Z];
}

// =====================================================================
// NuclearDensity class
// Implements the nucdens() function with persistent per-irho state.
// =====================================================================
class NuclearDensity {
public:
    static constexpr int MAX_IRHO = 10;

    NuclearDensity() {
        for (int i = 0; i < MAX_IRHO; ++i) {
            firsttime_[i] = true;
            Asave_[i] = 0;
            Zsave_[i] = 0;
            NWS_[i] = 0.0;
            RA_[i] = 0.0;
            RAequiv_[i] = 0.0;
            aa_[i] = 0.0;
            cc_[i] = 0.0;
        }
    }

    // Main density function.
    // r     = distance from nuclear center (fm)
    // Z     = atomic number
    // A     = atomic mass
    // idist = 0: Hard Sphere, 1: Woods-Saxon
    // irho  = density index (1-based, as in Fortran; mapped to 0-based internally)
    // Returns: nuclear density normalized to 1
    // Also updates RAeq, RAws output state.
    double operator()(double r, int Z, int A, int idist, int irho) {
        return compute(r, Z, A, idist, irho);
    }

    double compute(double r, int Z, int A, int idist, int irho) {
        static constexpr double onethird = 0.333333333333333;
        static constexpr double srft = 1.29099;
        static constexpr double threefourthoverpi = 0.2387336394417;

        int idx = irho - 1; // Convert to 0-based

        if (!firsttime_[idx]) {
            if (Asave_[idx] != A || Zsave_[idx] != Z) {
                firsttime_[idx] = true;
            }
        }
        Asave_[idx] = A;
        Zsave_[idx] = Z;

        double result;

        if (Z == 1 && A == 2) {
            // REID's SOFT-CORE (Deuterium)
            if (firsttime_[idx]) {
                RAequiv_[idx] = srft * std::sqrt(
                    dgauss11([](double rr) { return r2DD(rr); }, 0.0, 20.0, 1e-5));
                firsttime_[idx] = false;
            }
            result = DeuteronDensity(r);
            RAeq = RAequiv_[idx];
            RAws = 1.4111;

        } else if (idist == 0) {
            // HARD SPHERE
            RAeq = 1.12 * std::pow(static_cast<double>(A), onethird);
            RAws = RAeq;
            if (r <= RAeq) {
                result = threefourthoverpi / (RAeq * RAeq * RAeq);
            } else {
                result = 0.0;
            }

        } else if (idist == 1) {
            // WOODS-SAXON
            if (firsttime_[idx]) {
                WoodsSaxonParams ws_local;
                if (Z == 0 || Z > 92) {
                    // Bialas parametrization [4]
                    double Athird = std::pow(static_cast<double>(A), onethird);
                    ws_local.RR = (0.978 + 0.0206 * Athird) * Athird;
                    ws_local.a0 = 0.523;
                    ws_local.c0 = 0.0;
                } else {
                    WSParamEntry entry = get_ws_params(Z);
                    if (entry.R < 0.1) {
                        // No experimental data; use Bialas parametrization
                        double Athird = std::pow(static_cast<double>(A), onethird);
                        ws_local.RR = (0.978 + 0.0206 * Athird) * Athird;
                        ws_local.a0 = 0.523;
                        ws_local.c0 = 0.0;
                    } else {
                        ws_local.RR = entry.R;
                        ws_local.a0 = entry.a;
                        ws_local.c0 = entry.c;
                    }
                }
                ws_ = ws_local;

                double limit = 10.0 * ws_.RR;
                NWS_[idx] = 1.0 / dgauss11(
                    [this](double rr) { return WoodsSaxon3d(rr, ws_); },
                    0.0, limit, 1e-5);
                RAequiv_[idx] = srft * std::sqrt(NWS_[idx] *
                    dgauss11(
                        [this](double rr) { return r2WS3d(rr, ws_); },
                        0.0, limit, 1e-5));
                RA_[idx] = ws_.RR;
                aa_[idx] = ws_.a0;
                cc_[idx] = ws_.c0;
                firsttime_[idx] = false;
            }

            RAeq = RAequiv_[idx];
            RAws = RA_[idx];

            if (r - RA_[idx] < 700.0) {
                result = NWS_[idx]
                    * (1.0 + cc_[idx] * (r * r / (RA_[idx] * RA_[idx])))
                    / (1.0 + std::exp((r - RA_[idx]) / aa_[idx]));
            } else {
                result = 0.0;
            }
            if (result < 0.0) result = 0.0;

        } else {
            std::fprintf(stderr,
                "ERROR (nucdens): called out of range: idist,Z,A= %d %d %d\n",
                idist, Z, A);
            std::abort();
        }

        return result;
    }

    // Public output state (set after each call, mirrors Fortran common block)
    double RAeq = 0.0;   // hard-sphere equivalent radius
    double RAws = 0.0;   // Woods-Saxon radius parameter

private:
    bool firsttime_[MAX_IRHO];
    int Asave_[MAX_IRHO];
    int Zsave_[MAX_IRHO];
    double NWS_[MAX_IRHO];
    double RA_[MAX_IRHO];
    double RAequiv_[MAX_IRHO];
    double aa_[MAX_IRHO];
    double cc_[MAX_IRHO];
    WoodsSaxonParams ws_;  // current Woods-Saxon params for integration lambdas
};

// =====================================================================
// DensityTable: replaces GenNucDens subroutine
// Generates a cumulative density table for sampling interaction positions.
// =====================================================================
struct DensityTable {
    static constexpr int TABLE_SIZE = 2000;

    double density_table[TABLE_SIZE];
    double quantity_table[TABLE_SIZE];
    double step_size;
    double init_dens;

    DensityTable() : step_size(0.0), init_dens(0.0) {
        for (int i = 0; i < TABLE_SIZE; ++i) {
            density_table[i] = 0.0;
            quantity_table[i] = 0.0;
        }
    }

    // Generate the density table for a given nucleus.
    // iZ    = atomic number
    // iA    = atomic mass
    // iDens = distribution type (0=HS, 1=WS)
    void generate(int iZ, int iA, int iDens) {
        static constexpr double pi = 3.141592653;

        NuclearDensity nucdens;
        int idist = iDens;
        int irho = 1;
        step_size = 0.01;
        init_dens = 0.005;
        double r = init_dens;
        double integral = 0.0;

        for (int i = 0; i < TABLE_SIZE; ++i) {
            density_table[i] = nucdens.compute(r, iZ, iA, idist, irho);
            integral += 4.0 * pi * density_table[i] * r * r * step_size;
            quantity_table[i] = integral;
            r += step_size;
        }
    }
};

// =====================================================================
// InteractionPosition: replaces InterPos subroutine
// Samples a random interaction position from a density table.
// =====================================================================
struct InteractionPosition {
    double pos_radius;
    double pos_theta;
    double pos_phi;
    double x;
    double y;
    double z;

    InteractionPosition()
        : pos_radius(0.0), pos_theta(0.0), pos_phi(0.0),
          x(0.0), y(0.0), z(0.0) {}

    // Sample an interaction position from the density table using the
    // provided random number generator.
    // rng must provide operator() returning a double in [0,1).
    template <typename RNG>
    void sample(const DensityTable& table, RNG& rng) {
        static constexpr double pi = 3.141592653;

        // Random position in cumulative distribution
        pos_radius = rng() * table.quantity_table[DensityTable::TABLE_SIZE - 1];
        pos_theta = std::acos(2.0 * rng() - 1.0);
        pos_phi = 2.0 * pi * rng();

        // Find radius bin by walking the cumulative table
        double r = table.init_dens;
        int i = 0;
        while (i < DensityTable::TABLE_SIZE && pos_radius >= table.quantity_table[i]) {
            r += table.step_size;
            ++i;
        }
        pos_radius = r + (rng() - 0.5) * table.step_size;

        // Convert to Cartesian coordinates
        x = pos_radius * std::sin(pos_theta) * std::cos(pos_phi);
        y = pos_radius * std::sin(pos_theta) * std::sin(pos_phi);
        z = pos_radius * std::cos(pos_theta);
    }
};

} // namespace farm

#endif // FARM_NUCLEAR_DENSITY_H

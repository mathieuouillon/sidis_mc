#ifndef FARM_QUENCHING_H
#define FARM_QUENCHING_H

// ============================================================================
// QuenchingEngine: C++ translation of qweight.f90
//
// Implements Salgado-Wiedemann and Arleo quenching weight models.
// Original Fortran by C.A. Salgado, U.A. Wiedemann, A. Accardi, F. Arleo.
// ============================================================================

#include "pythia6.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace farm {

class QuenchingEngine {
public:
    // --- Configuration (mirrors quenching_module + config_module) -----------

    // Quenching weight model: 1 = Salgado-Wiedemann, 2 = Arleo
    int iqw   = 1;
    // Finite-size corrections flag (0 = no, 1 = yes)
    int scor  = 0;
    // Finite-energy corrections flag (0 = no, 1 = yes)
    int ncor  = 0;
    // 1 = multiple soft scatterings, 2 = single hard scattering
    int sfthrd = 1;

    // Strong coupling at soft scales
    double alphas = 1.0 / 3.0;

    // Config flags
    int iqg  = 0;   // include gluon quenching (1 = yes)
    int iEg  = 0;   // add energy-conservation gluon (1 = yes)
    int iPtF = 0;   // pT broadening model (0..3)

    // Transport coefficients
    double qhat   = 0.0;
    double ehat   = 0.0;
    double SupFac = 1.0;

    // --- Observables (output, mirrors QW_* variables) ----------------------

    int    QW_nb   = 0;
    double QW_w    = 0.0;
    double QW_L    = 0.0;
    double QW_wc   = 0.0;
    double QW_R    = 0.0;
    double QW_chi  = 0.0;
    double QW_th   = 0.0;
    double QW_qhat = 0.0;

    // --- Constructor -------------------------------------------------------

    QuenchingEngine() = default;

    // -----------------------------------------------------------------------
    // apply(): main entry point -- loops over PYTHIA particles, applies
    // energy loss.  Replaces ApplyQW().
    //
    //  py            -- Pythia6 event record
    //  x_inter, y_inter, z_inter -- interaction vertex position
    //  density_table -- nuclear density table (indexed by int(r/step))
    //  step_size_dens -- step size used to index density_table
    //  n_density     -- number of entries in density_table
    // -----------------------------------------------------------------------
    void apply(Pythia6& py,
               double x_inter, double y_inter, double z_inter,
               const double* density_table, double step_size_dens, int /*n_density*/)
    {
        const double cutoff = 0.4;

        int ip = 1;
        while (ip <= py.n_particles()) {
            int kid = py.get_k(ip, 2);
            int kst = py.get_k(ip, 1);
            double E = py.get_p(ip, 4);

            bool is_quark = (std::abs(kid) < 4);
            bool is_gluon = (kid == 21 && iqg == 1);

            if ((is_quark || is_gluon) && kst < 9 && E > cutoff) {

                // Stock initial values
                double inix = py.get_p(ip, 1);
                double iniy = py.get_p(ip, 2);
                double iniz = py.get_p(ip, 3);

                // Normalised initial direction
                double tot = std::sqrt(inix * inix + iniy * iniy + iniz * iniz);
                double ipix = inix / tot;
                double ipiy = iniy / tot;
                double ipiz = iniz / tot;

                // Compute quenching weight
                QW_nb += 1;
                qw_comput(py.get_p(ip, 1), py.get_p(ip, 2), py.get_p(ip, 3),
                          py.get_p(ip, 4), py.get_k(ip, 2),
                          x_inter, y_inter, z_inter,
                          density_table, step_size_dens, py);

                if (QW_w > 0.0) {
                    double ipt = 0.0;

                    // Determine transverse momentum of final parton
                    if (iPtF == 0) {
                        ipt = 0.0;
                    } else if (iPtF == 1) {
                        ipt = qhat * QW_L;
                    } else if (iPtF == 2) {
                        ipt = 4.0 * QW_w / 3.0 / alphas / QW_L * SupFac * SupFac;
                    } else if (iPtF == 3) {
                        ipt = (QW_w * std::sin(QW_th) * SupFac)
                            * (QW_w * std::sin(QW_th) * SupFac);
                    }

                    double ipl = 0.0;
                    E = py.get_p(ip, 4);

                    // Implement energy loss and pT
                    if (E - QW_w < cutoff) {
                        double th = py.random() * 2.0 * 3.14159265 - 3.14159265;
                        ipl = std::cos(th) * cutoff;
                        ipt = std::sin(th) * cutoff;
                    } else {
                        double iptot = (E - QW_w) * (E - QW_w);
                        if (iptot > ipt) {
                            ipl = iptot - ipt;
                            ipt = std::sqrt(ipt);
                            ipl = std::sqrt(ipl);
                        } else {
                            ipl = 0.0;
                            ipt = std::sqrt(iptot);
                        }
                    }

                    // Generate normalised transverse vector
                    double ph = 4.0 * std::asin(1.0) * py.random();
                    double iptx = (ipiz - ipiy) * std::cos(ph)
                                - (ipix * ipiy + ipix * ipiz) * std::sin(ph);
                    double ipty = ipix * std::cos(ph)
                                + (ipix * ipix + ipiz * ipiz - ipiy * ipiz) * std::sin(ph);
                    double iptz = -ipix * std::cos(ph)
                                + (ipix * ipix + ipiy * ipiy - ipiy * ipiz) * std::sin(ph);
                    double tnorm = std::sqrt(iptx * iptx + ipty * ipty + iptz * iptz);
                    iptx /= tnorm;
                    ipty /= tnorm;
                    iptz /= tnorm;

                    // Check perpendicularity
                    double sca = ipix * iptx + ipiy * ipty + ipiz * iptz;
                    if (std::abs(sca) > 1.0e-6) {
                        std::cerr << "problem tot = " << sca << "\n";
                    }

                    // Generate new parton momenta
                    double ipx = ipt * iptx + ipl * ipix;
                    double ipy = ipt * ipty + ipl * ipiy;
                    double ipz = ipt * iptz + ipl * ipiz;

                    // qhat summing
                    double p1 = py.get_p(ip, 1);
                    double p2 = py.get_p(ip, 2);
                    double p3 = py.get_p(ip, 3);
                    double dot = ipx * p1 + ipy * p2 + ipz * p3;
                    double mag2 = p1 * p1 + p2 * p2 + p3 * p3;
                    QW_qhat += ((ipx * ipx + ipy * ipy + ipz * ipz) - dot * dot / mag2) / QW_L;

                    // Fill Pythia array
                    py.set_p(ip, 1, ipx);
                    py.set_p(ip, 2, ipy);
                    py.set_p(ip, 3, ipz);
                    double mass = py.get_p(ip, 5);
                    py.set_p(ip, 4, std::sqrt(mass * mass + ipx * ipx + ipy * ipy + ipz * ipz));

                    // Add gluon if requested
                    if (iEg == 1) {
                        double ipg = tot - py.get_p(ip, 4);
                        double ptg, plg;
                        if (ipt < ipg) {
                            ptg = -ipt;
                            plg = std::sqrt(ipg * ipg - ptg * ptg);
                        } else {
                            ptg = ipg;
                            plg = 0.0;
                        }

                        double ipgx = ptg * iptx + plg * ipix;
                        double ipgy = ptg * ipty + plg * ipiy;
                        double ipgz = ptg * iptz + plg * ipiz;

                        // Shift particles to make room for the gluon
                        int N = py.n_particles();
                        for (int iq = ip; iq <= N; ++iq) {
                            int ir = ip + N - iq;
                            if (py.get_k(ip - 1, 1) == 2 || ir != ip) {
                                for (int j = 1; j <= 5; ++j) {
                                    py.set_p(ir + 1, j, py.get_p(ir, j));
                                    py.set_k(ir + 1, j, py.get_k(ir, j));
                                }
                            }
                        }

                        if (py.get_k(ip - 1, 1) == 2) {
                            py.set_p(ip, 1, ipgx);
                            py.set_p(ip, 2, ipgy);
                            py.set_p(ip, 3, ipgz);
                            py.set_p(ip, 4, ipg);
                            py.set_p(ip, 5, 0.0);

                            py.set_k(ip, 1, 2);
                            py.set_k(ip, 2, 21);
                            py.set_k(ip, 3, ip);
                            ip += 1;
                        } else {
                            ip += 1;
                            py.set_p(ip, 1, ipgx);
                            py.set_p(ip, 2, ipgy);
                            py.set_p(ip, 3, ipgz);
                            py.set_p(ip, 4, ipg);
                            py.set_p(ip, 5, 0.0);

                            py.set_k(ip, 1, 2);
                            py.set_k(ip, 2, 21);
                            py.set_k(ip, 3, ip);
                        }
                        py.set_n(N + 1);
                    }
                } // QW_w > 0
            } // parton selection
            ip += 1;
        } // while
    }

private:
    // --- Lookup tables for SW quenching weights ----------------------------

    // Multiple soft scattering tables
    double mult_xx[400]  = {};
    double mult_daq[34]  = {};
    double mult_caq[34][261] = {};
    double mult_rrr[34]  = {};
    double mult_xxg[400] = {};
    double mult_dag[34]  = {};
    double mult_cag[34][261] = {};
    double mult_rrrg[34] = {};

    // Single hard scattering tables
    double lin_xx[400]   = {};
    double lin_daq[34]   = {};
    double lin_caq[34][261] = {};
    double lin_rrr[34]   = {};
    double lin_xxg[400]  = {};
    double lin_dag[34]   = {};
    double lin_cag[34][261] = {};
    double lin_rrrg[34]  = {};

    // Initialization flags
    bool mult_initialized = false;
    bool lin_initialized  = false;

    // Reweighting flag (per-event)
    int irw = 0;

    // -----------------------------------------------------------------------
    // qw_comput(): compute quenching weight for one parton.
    // Replaces QWComput().
    // -----------------------------------------------------------------------
    void qw_comput(double ipx, double ipy, double ipz, double E, int id,
                   double x_inter, double y_inter, double z_inter,
                   const double* density_table, double step_size_dens,
                   Pythia6& py)
    {
        QW_w  = 0.0;
        QW_L  = 0.0;
        QW_wc = 0.0;
        QW_R  = 0.0;
        QW_th = 0.0;

        double d = 0.0;
        double qhateff = qhat + ehat;

        // Determine parton type: 0 = gluon, 1 = quark
        int ipart;
        if (id == 21) {
            ipart = 0;
        } else if (std::abs(id) < 7) {
            ipart = 1;
        } else {
            std::cerr << "Unknown parton with id = " << id << "\n";
            ipart = 1; // fallback
        }
        irw = 0;

        const int nb_step = 200;
        const double integral_step = 0.1;

        // Normalise momentum direction
        double pp = std::sqrt(ipx * ipx + ipy * ipy + ipz * ipz);
        double px = ipx / pp * integral_step;
        double py_step = ipy / pp * integral_step;
        double pz = ipz / pp * integral_step;

        double x = x_inter;
        double y = y_inter;
        double z = z_inter;

        // Integration to calculate wc and R
        double radius = std::sqrt(x * x + y * y + z * z);
        while (radius < 20.0) {
            int idx = static_cast<int>(radius / step_size_dens);
            QW_wc += integral_step * d * density_table[idx];
            QW_R  += integral_step * density_table[idx];
            d += integral_step;
            x += px;
            y += py_step;
            z += pz;
            radius = std::sqrt(x * x + y * y + z * z);
        }

        QW_L  = QW_wc / QW_R;
        QW_wc = qhateff / density_table[0] * QW_wc;
        QW_R  = 2.0 * density_table[0] * QW_wc * QW_wc / QW_R / qhateff;

        // Convert units: fm -> GeV^-1, GeV^2*fm -> GeV, GeV^2*fm^2 -> dimensionless
        QW_wc = QW_wc / 0.1973269;
        QW_R  = QW_R / (0.1973269 * 0.1973269);

        // Calculate energy loss probability
        double step_QW;
        if (sfthrd == 1) step_QW = 2.5 / nb_step;
        if (sfthrd == 2) step_QW = 9.8 / nb_step;

        double yy = E / QW_wc;

        double cont[1000];
        double disc = 0.0;
        double total = 0.0;

        for (int i = 1; i <= nb_step; ++i) {
            double xx = step_QW * i;
            qweight_calc(ipart, QW_R, xx, yy, cont[i - 1], disc);
            total += cont[i - 1] * step_QW;
        }
        total += disc;
        disc /= total;
        for (int i = 0; i < nb_step; ++i) {
            cont[i] /= total;
        }

        // Pick randomly a quenching weight from the table
        if (disc < 1.0) {
            double randnum = py.random();
            if (randnum > disc) {
                total = disc;
                int i = 0;
                while (randnum > total) {
                    total += cont[i] * step_QW;
                    i += 1;
                }
                QW_w = i * step_QW * QW_wc;
            }
        }

        // Calculate the angle probability
        if (QW_w > 0.0) {
            step_QW = 1.0 / nb_step;
            yy = E / QW_wc;
            double xx_w = QW_w / QW_wc;

            total = 0.0;
            for (int i = 1; i <= nb_step; ++i) {
                double ChiR = (step_QW * i) * (step_QW * i) * QW_R;
                qweight_calc(ipart, ChiR, xx_w, yy, cont[i - 1], disc);
                // Do not keep negative probabilities
                if (cont[i - 1] < 0.0) cont[i - 1] = 0.0;
                total += cont[i - 1] * step_QW;
            }
            for (int i = 0; i < nb_step; ++i) {
                cont[i] /= total;
            }
            double randnum = py.random();
            total = 0.0;
            int i = 0;
            while (randnum > total) {
                total += cont[i] * step_QW;
                i += 1;
            }
            QW_chi = i * step_QW;
            QW_th  = std::asin(QW_chi);
        }

        if (std::isnan(QW_th)) QW_th = 3.14159 / 2.0;
    }

    // -----------------------------------------------------------------------
    // qweight_calc(): dispatcher for quenching weight models.
    // Replaces the Fortran qweight() subroutine.
    // -----------------------------------------------------------------------
    void qweight_calc(int ipart, double rrrr, double xx, double yy,
                      double& cont, double& disc)
    {
        if (scor == 1 && ncor == 1) {
            throw std::runtime_error(
                "qweight: finite size & finite energy corrections not yet implemented");
        }

        // Arleo asymptotic medium size
        if (iqw == 2) {
            disc = 0.0;
            if (alphas <= 0.0) {
                throw std::runtime_error("qweight: alphas < 0");
            }
            double kk = 1.0 / (2.0 * alphas);
            if (ipart == 0) {
                cont = kk * dbarg(kk * xx, yy);
            } else {
                cont = kk * dbarq(kk * xx, yy);
            }

        // Salgado-Wiedemann quenching weights
        } else if (iqw == 1) {
            // Initialise tables on first call
            if (sfthrd == 1 && !mult_initialized) {
                initmult(alphas);
                mult_initialized = true;
            }
            if (sfthrd == 2 && !lin_initialized) {
                initlin(alphas);
                lin_initialized = true;
            }

            // If no size corrections, use very large R
            double rrin = (scor == 0) ? 1.0e8 : rrrr;

            static constexpr double xxmultmax = 2.59;
            static constexpr double xxlinmax  = 9.87;

            if (sfthrd == 1) {
                if (xx <= xxmultmax) {
                    swqmult(ipart, rrin, xx, cont, disc);
                } else {
                    swqmult(ipart, rrin, xxmultmax, cont, disc);
                    cont = 0.0;
                }
            } else {
                if (xx <= xxlinmax) {
                    swqlin(ipart, rrin, xx, cont, disc);
                } else {
                    swqlin(ipart, rrin, xxlinmax, cont, disc);
                    cont = 0.0;
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Arleo quenching weight functions
    // -----------------------------------------------------------------------

    double dbarg(double wl, double e) const {
        return (4.0 / 9.0) * dbarq((4.0 / 9.0) * wl, e);
    }

    double dbarq(double wl, double e) const {
        static constexpr double pi = 3.1415926;
        if (wl == 0.0) {
            return 0.0;
        }
        double mu  = xmu(e);
        double sig = xsigma(e);
        return std::exp(-(std::log(wl) - mu) * (std::log(wl) - mu)
                        / (2.0 * sig * sig))
               / (std::sqrt(2.0 * pi) * sig * wl);
    }

    double xmu(double e) const {
        if (ncor == 0) {
            return -1.5;
        } else {
            return -1.5 + 0.81 * (std::exp(-0.2 / e) - 1.0);
        }
    }

    double xsigma(double e) const {
        if (ncor == 0) {
            return 0.72;
        } else {
            return 0.72 + 0.33 * (std::exp(-0.2 / e) - 1.0);
        }
    }

    // -----------------------------------------------------------------------
    // swqmult(): Salgado-Wiedemann interpolation for multiple soft scattering
    // -----------------------------------------------------------------------
    void swqmult(int ipart, double rrrr, double xxxx,
                 double& continuous, double& discrete)
    {
        continuous = 0.0;
        discrete   = 0.0;

        double rrin = rrrr;
        double xxin = xxxx;

        int nrlow = 0, nrhigh = 0;
        double rrhigh = 0.0, rrlow = 0.0;

        for (int nr = 0; nr < 34; ++nr) {
            if (rrin < mult_rrr[nr]) {
                rrhigh = mult_rrr[nr];
            } else {
                rrhigh = mult_rrr[nr - 1];
                rrlow  = mult_rrr[nr];
                nrlow  = nr;
                nrhigh = nr - 1;
                break;
            }
        }

        double rfraclow  = (rrhigh - rrin) / (rrhigh - rrlow);
        double rfrachigh = (rrin - rrlow) / (rrhigh - rrlow);
        if (rrin > 10000.0) {
            rfraclow  = std::log(rrhigh / rrin) / std::log(rrhigh / rrlow);
            rfrachigh = std::log(rrin / rrlow)  / std::log(rrhigh / rrlow);
        }

        // Quark: ipart != 0
        if (ipart != 0 && rrin >= mult_rrr[0]) {
            nrlow    = 0;
            nrhigh   = 0;
            rfraclow = 1.0;
            rfrachigh = 0.0;
        }
        // Gluon: ipart == 0
        if (ipart == 0 && rrin >= mult_rrrg[0]) {
            nrlow    = 0;
            nrhigh   = 0;
            rfraclow = 1.0;
            rfrachigh = 0.0;
        }

        // mult_xx[259] corresponds to Fortran mult_xx(260)
        if (xxxx >= mult_xx[259]) {
            std::cerr << "swqmult: xxxx >= mult_xx(260), stopping.\n";
            std::abort();
        }

        int nxlow  = static_cast<int>(xxin / 0.01);  // 0-based index
        int nxhigh = nxlow + 1;
        double xfraclow  = (mult_xx[nxhigh] - xxin) / 0.01;
        double xfrachigh = (xxin - mult_xx[nxlow]) / 0.01;

        double clow, chigh;
        if (ipart != 0) {
            clow  = xfraclow * mult_caq[nrlow][nxlow]  + xfrachigh * mult_caq[nrlow][nxhigh];
            chigh = xfraclow * mult_caq[nrhigh][nxlow] + xfrachigh * mult_caq[nrhigh][nxhigh];
        } else {
            clow  = xfraclow * mult_cag[nrlow][nxlow]  + xfrachigh * mult_cag[nrlow][nxhigh];
            chigh = xfraclow * mult_cag[nrhigh][nxlow] + xfrachigh * mult_cag[nrhigh][nxhigh];
        }

        continuous = rfraclow * clow + rfrachigh * chigh;

        if (ipart != 0) {
            discrete = rfraclow * mult_daq[nrlow] + rfrachigh * mult_daq[nrhigh];
        } else {
            discrete = rfraclow * mult_dag[nrlow] + rfrachigh * mult_dag[nrhigh];
        }
    }

    // -----------------------------------------------------------------------
    // swqlin(): Salgado-Wiedemann interpolation for single hard scattering
    // -----------------------------------------------------------------------
    void swqlin(int ipart, double rrrr, double xxxx,
                double& continuous, double& discrete)
    {
        continuous = 0.0;
        discrete   = 0.0;

        double rrin = rrrr;
        double xxin = xxxx;

        int nrlow = 0, nrhigh = 0;
        double rrhigh = 0.0, rrlow = 0.0;

        for (int nr = 0; nr < 34; ++nr) {
            if (rrin < lin_rrr[nr]) {
                rrhigh = lin_rrr[nr];
            } else {
                rrhigh = lin_rrr[nr - 1];
                rrlow  = lin_rrr[nr];
                nrlow  = nr;
                nrhigh = nr - 1;
                break;
            }
        }

        double rfraclow  = (rrhigh - rrin) / (rrhigh - rrlow);
        double rfrachigh = (rrin - rrlow) / (rrhigh - rrlow);
        if (rrin > 10000.0) {
            rfraclow  = std::log(rrhigh / rrin) / std::log(rrhigh / rrlow);
            rfrachigh = std::log(rrin / rrlow)  / std::log(rrhigh / rrlow);
        }

        // Note: swqlin uses ipart==1 for quark (unlike swqmult which uses ipart!=0)
        if (ipart == 1 && rrin >= lin_rrr[0]) {
            nrlow    = 0;
            nrhigh   = 0;
            rfraclow = 1.0;
            rfrachigh = 0.0;
        }
        if (ipart != 1 && rrin >= lin_rrrg[0]) {
            nrlow    = 0;
            nrhigh   = 0;
            rfraclow = 1.0;
            rfrachigh = 0.0;
        }

        // lin_xx[259] corresponds to Fortran lin_xx(260)
        if (xxxx < lin_xx[259]) {
            int nxlow_l  = static_cast<int>(xxin / 0.038);
            int nxhigh_l = nxlow_l + 1;
            double xfraclow  = (lin_xx[nxhigh_l] - xxin) / 0.038;
            double xfrachigh = (xxin - lin_xx[nxlow_l]) / 0.038;

            double clow, chigh;
            if (ipart == 1) {
                clow  = xfraclow * lin_caq[nrlow][nxlow_l]  + xfrachigh * lin_caq[nrlow][nxhigh_l];
                chigh = xfraclow * lin_caq[nrhigh][nxlow_l] + xfrachigh * lin_caq[nrhigh][nxhigh_l];
            } else {
                clow  = xfraclow * lin_cag[nrlow][nxlow_l]  + xfrachigh * lin_cag[nrlow][nxhigh_l];
                chigh = xfraclow * lin_cag[nrhigh][nxlow_l] + xfrachigh * lin_cag[nrhigh][nxhigh_l];
            }

            continuous = rfraclow * clow + rfrachigh * chigh;
        }

        if (ipart == 1) {
            discrete = rfraclow * lin_daq[nrlow] + rfrachigh * lin_daq[nrhigh];
        } else {
            discrete = rfraclow * lin_dag[nrlow] + rfrachigh * lin_dag[nrhigh];
        }
    }

    // -----------------------------------------------------------------------
    // initmult(): load lookup tables for multiple soft scattering
    //
    // Files: datafiles/qweight/cont0{3,5}.all, disc0{3,5}.all
    // Each cont line: xx  caq(1..34)  ->  35 columns
    // Each disc line: rrr  daq        ->  2 columns
    // First 261 lines = quark, next 261 lines = gluon (continuous)
    // First 34 lines = quark, next 34 lines = gluon (discrete)
    // -----------------------------------------------------------------------
    void initmult(double as) {
        std::string cont_file, disc_file;

        if (std::lround(as * 3.0) == 1) {
            cont_file = "datafiles/qweight/cont03.all";
            disc_file = "datafiles/qweight/disc03.all";
        } else if (std::lround(as * 2.0) == 1) {
            cont_file = "datafiles/qweight/cont05.all";
            disc_file = "datafiles/qweight/disc05.all";
        } else {
            throw std::runtime_error("initmult: alphas must be 1/3 or 1/2");
        }

        // Read continuous weights
        {
            std::ifstream fin(cont_file);
            if (!fin.is_open()) {
                throw std::runtime_error("Cannot open " + cont_file);
            }
            // Quark: 261 lines, each with xx + 34 values
            for (int nn = 0; nn < 261; ++nn) {
                fin >> mult_xx[nn];
                for (int j = 0; j < 34; ++j) {
                    fin >> mult_caq[j][nn];
                }
            }
            // Gluon: 261 lines
            for (int nn = 0; nn < 261; ++nn) {
                fin >> mult_xxg[nn];
                for (int j = 0; j < 34; ++j) {
                    fin >> mult_cag[j][nn];
                }
            }
        }

        // Read discrete weights
        {
            std::ifstream fin(disc_file);
            if (!fin.is_open()) {
                throw std::runtime_error("Cannot open " + disc_file);
            }
            // Quark: 34 lines
            for (int nn = 0; nn < 34; ++nn) {
                fin >> mult_rrr[nn] >> mult_daq[nn];
            }
            // Gluon: 34 lines
            for (int nn = 0; nn < 34; ++nn) {
                fin >> mult_rrrg[nn] >> mult_dag[nn];
            }
        }
    }

    // -----------------------------------------------------------------------
    // initlin(): load lookup tables for single hard scattering
    //
    // Files: datafiles/qweight/contlin0{3}.all, disclin0{3}.all
    // Same format as initmult.
    // -----------------------------------------------------------------------
    void initlin(double as) {
        std::string cont_file, disc_file;

        if (std::lround(as * 3.0) == 1) {
            cont_file = "datafiles/qweight/contlin03.all";
            disc_file = "datafiles/qweight/disclin03.all";
        } else if (std::lround(as * 2.0) == 1) {
            throw std::runtime_error("initlin: alphas=0.5 not yet implemented");
        } else {
            throw std::runtime_error("initlin: alphas must be 1/3 or 1/2");
        }

        // Read continuous weights
        {
            std::ifstream fin(cont_file);
            if (!fin.is_open()) {
                throw std::runtime_error("Cannot open " + cont_file);
            }
            // Quark: 261 lines
            for (int nn = 0; nn < 261; ++nn) {
                fin >> lin_xx[nn];
                for (int j = 0; j < 34; ++j) {
                    fin >> lin_caq[j][nn];
                }
            }
            // Gluon: 261 lines
            for (int nn = 0; nn < 261; ++nn) {
                fin >> lin_xxg[nn];
                for (int j = 0; j < 34; ++j) {
                    fin >> lin_cag[j][nn];
                }
            }
        }

        // Read discrete weights
        {
            std::ifstream fin(disc_file);
            if (!fin.is_open()) {
                throw std::runtime_error("Cannot open " + disc_file);
            }
            // Quark: 34 lines
            for (int nn = 0; nn < 34; ++nn) {
                fin >> lin_rrr[nn] >> lin_daq[nn];
            }
            // Gluon: 34 lines
            for (int nn = 0; nn < 34; ++nn) {
                fin >> lin_rrrg[nn] >> lin_dag[nn];
            }
        }
    }
};

} // namespace farm

#endif // FARM_QUENCHING_H

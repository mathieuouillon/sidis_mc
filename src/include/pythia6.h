#ifndef FARM_PYTHIA6_H
#define FARM_PYTHIA6_H

// ============================================================================
// Pythia6: C++ wrapper around PYTHIA 6.4.28 Fortran library
//
// Accesses PYTHIA routines and common blocks via extern "C" linkage.
// Fortran common blocks are memory-mapped as C structs.
// Array indexing: Fortran P(i,j) = C P[j-1][i-1] (column-major match).
// ============================================================================

#include <cstring>

// --- Fortran PYTHIA6 routine declarations -----------------------------------
extern "C" {
    void pyinit_(const char* frame, const char* beam, const char* target,
                 double* energy, int frame_len, int beam_len, int target_len);
    void pyevnt_();
    void pyexec_();
    double pyr_(int* dummy);

    // Common blocks (Fortran column-major layout)
    extern struct { int    MSTP[200]; double PARP[200]; int MSTI[200]; double PARI[200]; } pypars_;
    extern struct { int    MSTU[200]; double PARU[200]; int MSTJ[200]; double PARJ[200]; } pydat1_;
    extern struct { int    N, NPAD, K[5][4000]; double P[5][4000], V[5][4000]; } pyjets_;
    extern struct { int    MSEL, MSELPD, MSUB[500]; int KFIN[2][81]; double CKIN[200]; } pysubs_;
    extern struct { int    NGENPD; int NGEN[3][501]; double XSEC[3][501]; } pyint5_;
    extern struct { int    MRPY[6]; double RRPY[100]; } pydatr_;
}

// --- Pythia6 C++ class ------------------------------------------------------

class Pythia6 {
public:
    // --- Initialization -----------------------------------------------------

    void init_proton(double energy) {
        pyinit_("FIXT", "gamma/e-", "p+", &energy, 4, 8, 2);
    }

    void init_neutron(double energy) {
        pyinit_("FIXT", "gamma/e-", "n0", &energy, 4, 8, 2);
    }

    // --- Event generation ---------------------------------------------------

    void generate_event() { pyevnt_(); }
    void fragment()       { pyexec_(); }

    double random() {
        int dummy = 0;
        return pyr_(&dummy);
    }

    // --- Configuration (PYPARS) ---------------------------------------------

    void set_mstp(int i, int val)    { pypars_.MSTP[i-1] = val; }
    int  get_mstp(int i) const       { return pypars_.MSTP[i-1]; }
    void set_parp(int i, double val) { pypars_.PARP[i-1] = val; }
    double get_parp(int i) const     { return pypars_.PARP[i-1]; }
    int  get_msti(int i) const       { return pypars_.MSTI[i-1]; }
    double get_pari(int i) const     { return pypars_.PARI[i-1]; }

    // --- Configuration (PYDAT1) ---------------------------------------------

    void set_mstj(int i, int val)    { pydat1_.MSTJ[i-1] = val; }
    int  get_mstj(int i) const       { return pydat1_.MSTJ[i-1]; }
    void set_parj(int i, double val) { pydat1_.PARJ[i-1] = val; }
    double get_parj(int i) const     { return pydat1_.PARJ[i-1]; }
    void set_mstu(int i, int val)    { pydat1_.MSTU[i-1] = val; }
    int  get_mstu(int i) const       { return pydat1_.MSTU[i-1]; }
    void set_paru(int i, double val) { pydat1_.PARU[i-1] = val; }
    double get_paru(int i) const     { return pydat1_.PARU[i-1]; }

    // --- Configuration (PYSUBS) ---------------------------------------------

    void set_msel(int val)           { pysubs_.MSEL = val; }
    int  get_msel() const            { return pysubs_.MSEL; }
    void set_ckin(int i, double val) { pysubs_.CKIN[i-1] = val; }
    double get_ckin(int i) const     { return pysubs_.CKIN[i-1]; }

    // --- Fragmentation control ----------------------------------------------

    void disable_fragmentation() { set_mstj(1, 0); }
    void enable_fragmentation()  { set_mstj(1, 1); }

    // --- Event record (PYJETS) — 1-indexed like Fortran ---------------------

    int n_particles() const { return pyjets_.N; }
    void set_n(int n)       { pyjets_.N = n; }

    // K(i,j): status/flavor array
    int  get_k(int i, int j) const     { return pyjets_.K[j-1][i-1]; }
    void set_k(int i, int j, int val)  { pyjets_.K[j-1][i-1] = val; }
    int  status(int i) const           { return get_k(i, 1); }
    int  pdg_id(int i) const           { return get_k(i, 2); }
    int  mother(int i) const           { return get_k(i, 3); }

    // P(i,j): 4-momentum + mass
    double get_p(int i, int j) const      { return pyjets_.P[j-1][i-1]; }
    void   set_p(int i, int j, double val){ pyjets_.P[j-1][i-1] = val; }
    double px(int i) const     { return get_p(i, 1); }
    double py(int i) const     { return get_p(i, 2); }
    double pz(int i) const     { return get_p(i, 3); }
    double energy(int i) const { return get_p(i, 4); }
    double mass(int i) const   { return get_p(i, 5); }

    // Add a particle to the event record, returns the new particle index
    int add_particle(int stat, int id, int moth,
                      double ppx, double ppy, double ppz, double e, double m) {
        int n = pyjets_.N + 1;
        pyjets_.N = n;
        set_k(n, 1, stat);
        set_k(n, 2, id);
        set_k(n, 3, moth);
        set_p(n, 1, ppx);
        set_p(n, 2, ppy);
        set_p(n, 3, ppz);
        set_p(n, 4, e);
        set_p(n, 5, m);
        return n;
    }

    // --- Cross section (PYINT5) ---------------------------------------------

    double xsec(int isub, int j = 1) const { return pyint5_.XSEC[j-1][isub]; }

    // --- RNG seeding (PYDATR) -----------------------------------------------

    void seed(int iseed) {
        pydatr_.MRPY[1] = 0;           // MRPY(2)
        pydatr_.MRPY[2] = iseed % 85635; // MRPY(3)
        pydatr_.MRPY[3] = iseed % 67;    // MRPY(4)
        pydatr_.MRPY[4] = iseed % 56;    // MRPY(5)
    }

    // --- DIS configuration (replaces PythiaConfigDIS in pythiaconfig.f90) ----

    void configure_dis() {
        set_mstp(14, 26);    // DIS only (photon structure)
        set_mstp(82, 0);     // No multiple interactions
        set_parp(2, 2.0);    // Min CM energy = 2 GeV
        set_parj(33, 0.3);   // Stop fragmentation threshold
        set_ckin(1, 2.0);    // Lower limit sqrt(s)
        set_ckin(3, 0.0);    // pT limits
        set_ckin(65, 1.0);   // Min Q2
        set_ckin(66, -1.0);  // Max Q2 (no limit)
        set_ckin(77, 2.0);   // Min W
        set_ckin(78, -1.0);  // Max W (no limit)
    }
};

#endif // FARM_PYTHIA6_H

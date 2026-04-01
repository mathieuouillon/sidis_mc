#ifndef FARM_INTERFACE_H
#define FARM_INTERFACE_H

// ============================================================================
// C interface to FARM Fortran physics routines (ISO_C_BINDING)
// ============================================================================

#ifdef __cplusplus
extern "C" {
#endif

// --- Configuration ----------------------------------------------------------
void farm_set_config(int nevent, int iTg, int iFM, float E0, int nkin, int seed);
void farm_set_physics_params(float FMlimit, int iIso, int iNS, int iLund,
                              int iQuenching, int iSim, int iqw,
                              float qhat, float ehat, int iDens,
                              int iqg, int iEg, int iPtF);
void farm_set_nucleon(int nucleon);
void farm_set_nkin(int nkin);
void farm_set_ievent(int ievent);
void farm_set_mstj1(int val);
void farm_set_vz(float vz);
void farm_set_nuclear_params(int iZ, int iA, float rFM);
void farm_init_physics(void);

// --- State access -----------------------------------------------------------
void farm_get_state(int* ievent, float* PPe, float* rFM,
                     int* iZ, int* iA, int* nucleon);
void farm_get_xsec99(double* xsec);
void farm_get_qw_stats(float* qhat, int* nb);
void farm_get_event_data(int* nb_part, int* ids, int* charges, int* mothers,
                          float* px, float* py, float* pz,
                          float* E, float* m, int max_part);

// --- Physics initialization -------------------------------------------------
void farm_init_nucl(void);
void farm_init_fm(void);
void farm_gen_nuc_dens(void);
void farm_init_random(void);
void farm_initialize(float* rFM, int* iZ, int* iA);

// --- Event generation -------------------------------------------------------
void farm_fm_param(void);
void farm_init_kin(void);
void farm_lorentz_fm(int transform_type);
void farm_lorentz_fm_back(int transform_type);
void farm_pythia_config_dis(void);
void farm_init_kin2book(void);
void farm_compute_v(void);
void farm_inter_pos(void);
void farm_apply_qw(void);
void farm_create_spec(void);
void farm_calc_vz(int target_type, float* vz);

// --- PYTHIA -----------------------------------------------------------------
void farm_pyinit_proton(double energy);
void farm_pyinit_neutron(double energy);
void farm_pyevnt(void);
void farm_pyexec(void);

// --- High-level event loop --------------------------------------------------
void farm_setup_kinematics(int ievent, int nevent);
void farm_generate_event(void);
void farm_needs_reinit(int ievent, int nkin_counter, int nevent, int* result);

// --- I/O --------------------------------------------------------------------
void farm_open_lund(const char* filename, int len, int* unit, int* err);
void farm_close_lund(int unit);
void farm_write_event(int unit);

// --- Timing -----------------------------------------------------------------
void farm_timex(double* t);

#ifdef __cplusplus
}
#endif

#endif // FARM_INTERFACE_H

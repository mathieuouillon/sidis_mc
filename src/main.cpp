#include <cstdio>
#include <cstdlib>
#include <string>
#include <random>
#include <chrono>
#include "include/farm_interface.h"
#include "include/physics.h"
#include "include/pythia6.h"
#include "include/event_processor.h"

// ============================================================================
// Configuration
// ============================================================================

struct SimConfig {
    // Required parameters
    int    nevent = 0;
    std::string output;
    int    target = -1;
    int    nkin   = 0;
    float  e0     = 0.0f;

    // Optional parameters
    int    iFM    = 5;
    int    seed   = -1;

    // Physics defaults
    float  FMlimit     = 0.5f;
    int    iIso        = 0;       // 0 = isospin symmetric, 1 = split
    int    iNS         = 0;       // nuclear spectators
    int    iQuenching  = 1;       // quenching weights
    int    iSim        = 1;       // PYTHIA simulation
    int    iqw         = 1;       // 1 = SW, 2 = Arleo
    float  qhat        = 0.34f;   // transport coefficient (GeV^2/fm)
    float  ehat        = 0.0f;    // drag coefficient
    int    iDens       = 1;       // 1 = Wood-Saxon
    int    iqg         = 1;       // quark and gluon quenching
    int    iEg         = 0;       // energy conservation gluon
    int    iPtF        = 3;       // pT broadening model

    static const char* target_name(int iTg) {
        static const char* names[] = {
            "p", "2H", "3H", "He3", "He4", "Li6", "Li7", "C",
            "Al", "Fe", "Sn", "Pb", "Ne", "Kr", "Xe", "Cu"
        };
        return (iTg >= 0 && iTg <= 15) ? names[iTg] : "Unknown";
    }
};

// ============================================================================
// CLI parsing
// ============================================================================

static SimConfig parse_args(int argc, char* argv[]) {
    SimConfig cfg;

    if (argc < 2) {
        printf("Usage: simulation --nevent <N> --output <file> --target <0-15> --nkin <N> --e0 <GeV>\n");
        printf("Optional: --iFM <0-5> --seed <N> --help\n");
        exit(1);
    }

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "--nevent" || arg == "-n") && i + 1 < argc)
            cfg.nevent = atoi(argv[++i]);
        else if ((arg == "--output" || arg == "-o") && i + 1 < argc)
            cfg.output = argv[++i];
        else if ((arg == "--target" || arg == "-t") && i + 1 < argc)
            cfg.target = atoi(argv[++i]);
        else if (arg == "--nkin" && i + 1 < argc)
            cfg.nkin = atoi(argv[++i]);
        else if (arg == "--e0" && i + 1 < argc)
            cfg.e0 = static_cast<float>(atof(argv[++i]));
        else if (arg == "--iFM" && i + 1 < argc)
            cfg.iFM = atoi(argv[++i]);
        else if ((arg == "--seed" || arg == "-s") && i + 1 < argc)
            cfg.seed = atoi(argv[++i]);
        else if (arg == "--help" || arg == "-h") {
            printf("Usage: simulation [OPTIONS]\n\n");
            printf("Required: --nevent <N>  --output <file>  --target <0-15>  --nkin <N>  --e0 <GeV>\n");
            printf("Optional: --iFM <0-5> [default: 5]  --seed <N> [default: time-based]\n\n");
            printf("Targets: 0=p 1=2H 2=3H 3=He3 4=He4 5=Li6 6=Li7 7=C\n");
            printf("         8=Al 9=Fe 10=Sn 11=Pb 12=Ne 13=Kr 14=Xe 15=Cu\n\n");
            printf("FM: 0=none 1=Bodek-Ritchie 2=Accardi-SVG 3=Accardi-CS 4=hard-sphere 5=Wiringa\n");
            exit(0);
        } else {
            fprintf(stderr, "Error: Unknown argument: %s\n", argv[i]);
            exit(1);
        }
    }

    if (cfg.nevent <= 0)    { fprintf(stderr, "Error: --nevent required\n"); exit(1); }
    if (cfg.output.empty()) { fprintf(stderr, "Error: --output required\n"); exit(1); }
    if (cfg.target < 0 || cfg.target > 15) { fprintf(stderr, "Error: --target 0-15 required\n"); exit(1); }
    if (cfg.nkin <= 0)      { fprintf(stderr, "Error: --nkin required\n"); exit(1); }
    if (cfg.e0 <= 0.0f)     { fprintf(stderr, "Error: --e0 required\n"); exit(1); }

    return cfg;
}

// ============================================================================
// LUND output writer (C++)
// ============================================================================

class LundWriter {
    FILE* file_ = nullptr;
public:
    bool open(const std::string& filename) {
        file_ = fopen(filename.c_str(), "w");
        return file_ != nullptr;
    }

    void write_event(const farm::EventData& evt, int iA, int iZ, float E0,
                      int nucleon, float vz) {
        if (!file_) return;
        int nb = static_cast<int>(evt.particles.size());
        fprintf(file_, "%12d%12d%12d%12d%12d%12d%12.7f%17d%12d%12.8f\n",
                nb, iA, iZ, 0, 0, 11, E0, nucleon, 1, 1.0f);
        for (int l = 0; l < nb; l++) {
            const auto& p = evt.particles[l];
            fprintf(file_, "%12d%12d%12d%12d%12d%12d%16.8E%16.8E%16.8E%16.8E%16.8E%15d%12d%12.8f\n",
                    l + 1, p.charge, 1, p.id, p.mother_id, 0,
                    p.px, p.py, p.pz, p.E, p.m, 0, 0, vz);
        }
    }

    void close() {
        if (file_) { fclose(file_); file_ = nullptr; }
    }

    ~LundWriter() { close(); }
};

// ============================================================================
// Simulation engine
// ============================================================================

class Simulation {
    SimConfig cfg_;
    Pythia6 pythia_;
    int iZ_ = 0, iA_ = 0;
    float rFM_ = 0.0f;
    int nkin_ = 0;
    std::mt19937 rng_;
    std::uniform_real_distribution<float> uniform_{0.0f, 1.0f};

public:
    explicit Simulation(const SimConfig& cfg)
        : cfg_(cfg), nkin_(cfg.nkin),
          rng_(cfg.seed >= 0 ? static_cast<unsigned>(cfg.seed) : std::random_device{}()) {}

    void initialize() {
        // Push configuration to Fortran modules
        farm_set_config(cfg_.nevent, cfg_.target, cfg_.iFM, cfg_.e0, cfg_.nkin, cfg_.seed);
        farm_set_physics_params(cfg_.FMlimit, cfg_.iIso, cfg_.iNS, 1 /*iLund*/,
                                cfg_.iQuenching, cfg_.iSim, cfg_.iqw,
                                cfg_.qhat, cfg_.ehat, cfg_.iDens,
                                cfg_.iqg, cfg_.iEg, cfg_.iPtF);

        // Skip frequent PYTHIA re-init for targets without Fermi motion
        if (cfg_.target == 0 || cfg_.iFM == 0) {
            nkin_ = cfg_.nevent + 1;
            farm_set_nkin(nkin_);
        }

        // Nuclear parameters lookup (C++)
        auto nuc = farm::get_nuclear_params(cfg_.target);
        iZ_ = nuc.Z;
        iA_ = nuc.A;
        rFM_ = nuc.rFM;

        // Push to Fortran modules and initialize FM tables, density, RNG
        farm_set_ievent(0);
        farm_set_nuclear_params(iZ_, iA_, rFM_);
        farm_init_physics();
    }

    void run(LundWriter& writer) {
        int nkin_counter = nkin_;

        for (int ievent = 0; ievent < cfg_.nevent; ievent++) {
            farm_set_ievent(ievent);

            // Re-initialize kinematics + PYTHIA when needed
            bool reinit = (ievent == 0)
                || (nkin_counter == nkin_)
                || (ievent == cfg_.nevent * iZ_ / iA_)
                || (ievent % 500000 == 0);

            if (reinit) {
                nkin_counter = 0;
                init_event_kinematics(ievent);
            }

            if (ievent % 10000 == 0)
                printf(" %d events processed\n", ievent);

            // Generate and process event
            auto evt = process_event();

            // Write LUND output (C++)
            float vz = farm::vertex_z(cfg_.target, uniform_(rng_), uniform_(rng_));
            float d1, d2, d3, d4;
            int nucleon_val;
            farm_get_fm_state(&d1, &d2, &d3, &d4, &nucleon_val);
            writer.write_event(evt, iA_, iZ_, cfg_.e0, nucleon_val, vz);

            nkin_counter++;
        }
    }

    void print_stats() const {
        printf(" X sec 99 =    %.16E\n", pythia_.xsec(99));

        float qw_qhat;
        int qw_nb;
        farm_get_qw_stats(&qw_qhat, &qw_nb);
        if (qw_nb > 0) printf(" q hat =    %E\n", qw_qhat / qw_nb);
    }

private:
    void init_event_kinematics(int ievent) {
        // Fermi motion + Lorentz boost (Fortran), returns beam energy
        double beam_energy;
        farm_setup_kinematics(ievent, cfg_.nevent, &beam_energy);

        // PYTHIA configuration + initialization (C++)
        pythia_.configure_dis();
        if (cfg_.iQuenching != 0)
            pythia_.disable_fragmentation();

        bool use_proton = (cfg_.iIso == 0)
            ? (ievent < cfg_.nevent * iZ_ / iA_)
            : (ievent < cfg_.nevent / 2);

        if (use_proton) {
            pythia_.init_proton(beam_energy);
        } else {
            printf(" X sec 99 =    %.16E\n", pythia_.xsec(99));
            pythia_.init_neutron(beam_energy);
        }
    }

    farm::EventData process_event() {
        // 1. PYTHIA event generation (C++)
        pythia_.generate_event();

        // 2. Inverse Lorentz boost back to lab frame (C++)
        if (cfg_.iFM != 0) {
            float BB1, B1x, B1y, B1z, Thi, Phi;
            farm_get_boost_params(&BB1, &B1x, &B1y, &B1z, &Thi, &Phi);
            farm::boost_all_back(pythia_, BB1, B1x, B1y, B1z, Thi, Phi);
        }

        // 3. Quenching (stays Fortran — complex physics)
        farm_apply_quenching();

        // 4. Fragmentation (C++ PYTHIA control)
        if (cfg_.iQuenching != 0) {
            pythia_.enable_fragmentation();
            pythia_.fragment();
            pythia_.disable_fragmentation();
        }

        // 5. Spectators (C++)
        if (cfg_.iNS == 1 && (cfg_.target >= 1 || cfg_.target <= 4)) {
            float nuc_the, nuc_phi, nuc_mom, FMintact;
            int nucleon;
            farm_get_fm_state(&nuc_the, &nuc_phi, &nuc_mom, &FMintact, &nucleon);
            farm::add_spectator(pythia_, cfg_.target, nucleon,
                                nuc_the, nuc_phi, nuc_mom, FMintact,
                                uniform_(rng_));
        }

        // 6. Compute DIS variables + extract particles (C++)
        return farm::compute_event(pythia_);
    }
};

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    setvbuf(stdout, NULL, _IOLBF, 0);

    SimConfig cfg = parse_args(argc, argv);

    printf(" Monte Carlo Simulation Parameters:\n");
    printf("   Number of events:  %d\n", cfg.nevent);
    printf("   Lund output file:  %s\n", cfg.output.c_str());
    printf("   Target type:       %s\n", SimConfig::target_name(cfg.target));
    printf("   Nkin value:        %d\n", cfg.nkin);
    printf("   Electron energy:   %.7f GeV\n", cfg.e0);
    printf("   Fermi motion:      %d\n", cfg.iFM);
    if (cfg.seed >= 0)
        printf("   Random seed:       %d\n", cfg.seed);
    else
        printf("   Random seed:       time-based\n");
    printf("\n");
    fflush(stdout);

    auto t1 = std::chrono::steady_clock::now();

    // Initialize simulation
    Simulation sim(cfg);
    sim.initialize();

    // Open output
    LundWriter writer;
    if (!writer.open(cfg.output)) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", cfg.output.c_str());
        return 1;
    }

    // Run event loop
    sim.run(writer);
    writer.close();

    // Statistics
    sim.print_stats();
    auto t2 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t2 - t1).count();
    printf(" %d events in %.2f s\n", cfg.nevent, elapsed);

    return 0;
}

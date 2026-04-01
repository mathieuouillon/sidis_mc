#include <cstdio>
#include <cstdlib>
#include <string>
#include <random>
#include "include/farm_interface.h"
#include "include/physics.h"

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
// LUND output writer (delegates to Fortran for format compatibility)
// ============================================================================

class LundWriter {
    int unit_ = -1;
    bool open_ = false;
public:
    bool open(const std::string& filename) {
        int err = 0;
        farm_open_lund(filename.c_str(), static_cast<int>(filename.size()), &unit_, &err);
        open_ = (err == 0);
        return open_;
    }

    void write_event() {
        if (open_) farm_write_event(unit_);
    }

    void close() {
        if (open_) { farm_close_lund(unit_); open_ = false; }
    }

    ~LundWriter() { close(); }
};

// ============================================================================
// Simulation engine
// ============================================================================

class Simulation {
    SimConfig cfg_;
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

        // Initialize nuclear physics, Fermi motion tables, density, and RNG
        farm_set_ievent(0);
        farm_initialize(&rFM_, &iZ_, &iA_);
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
                farm_setup_kinematics(ievent, cfg_.nevent);
            }

            // Progress
            if (ievent % 10000 == 0)
                printf(" %d events processed\n", ievent);

            // Generate and process one complete physics event
            farm_generate_event();

            // Vertex z position (C++) + write LUND output (Fortran I/O)
            float vz = farm::vertex_z(cfg_.target, uniform_(rng_), uniform_(rng_));
            farm_set_vz(vz);
            writer.write_event();

            nkin_counter++;
        }
    }

    void print_stats() const {
        double xsec;
        farm_get_xsec99(&xsec);
        printf(" X sec 99 =    %.16E\n", xsec);

        float qw_qhat;
        int qw_nb;
        farm_get_qw_stats(&qw_qhat, &qw_nb);
        if (qw_nb > 0) printf(" q hat =    %E\n", qw_qhat / qw_nb);
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

    double t1, t2;
    farm_timex(&t1);

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
    farm_timex(&t2);
    printf(" %d events in    %.16E s\n", cfg.nevent, t2 - t1);

    return 0;
}

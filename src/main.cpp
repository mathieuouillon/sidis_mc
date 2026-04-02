#include <array>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <random>
#include "include/physics.h"
#include "include/pythia6.h"
#include "include/event_processor.h"
#include "include/fermi_motion.h"
#include "include/accardi_fm.h"
#include "include/nuclear_density.h"
#include "include/quenching.h"

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

    static constexpr std::array<std::string_view, 16> TARGET_NAMES{
        "p", "2H", "3H", "He3", "He4", "Li6", "Li7", "C",
        "Al", "Fe", "Sn", "Pb", "Ne", "Kr", "Xe", "Cu"
    };
    static constexpr std::string_view target_name(int iTg) {
        return (iTg >= 0 && iTg < 16) ? TARGET_NAMES[iTg] : "Unknown";
    }
};

// ============================================================================
// CLI parsing
// ============================================================================

static SimConfig parse_args(int argc, char* argv[]) {
    SimConfig cfg;

    if (argc < 2) {
        throw std::invalid_argument(
            "Usage: simulation --nevent <N> --output <file> --target <0-15> --nkin <N> --e0 <GeV>\n"
            "Optional: --iFM <0-5> --seed <N> --help");
    }

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "--nevent" || arg == "-n") && i + 1 < argc)
            cfg.nevent = std::stoi(argv[++i]);
        else if ((arg == "--output" || arg == "-o") && i + 1 < argc)
            cfg.output = argv[++i];
        else if ((arg == "--target" || arg == "-t") && i + 1 < argc)
            cfg.target = std::stoi(argv[++i]);
        else if (arg == "--nkin" && i + 1 < argc)
            cfg.nkin = std::stoi(argv[++i]);
        else if (arg == "--e0" && i + 1 < argc)
            cfg.e0 = std::stof(argv[++i]);
        else if (arg == "--iFM" && i + 1 < argc)
            cfg.iFM = std::stoi(argv[++i]);
        else if ((arg == "--seed" || arg == "-s") && i + 1 < argc)
            cfg.seed = std::stoi(argv[++i]);
        else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: simulation [OPTIONS]\n\n"
                      << "Required: --nevent <N>  --output <file>  --target <0-15>  --nkin <N>  --e0 <GeV>\n"
                      << "Optional: --iFM <0-5> [default: 5]  --seed <N> [default: time-based]\n\n"
                      << "Targets: 0=p 1=2H 2=3H 3=He3 4=He4 5=Li6 6=Li7 7=C\n"
                      << "         8=Al 9=Fe 10=Sn 11=Pb 12=Ne 13=Kr 14=Xe 15=Cu\n\n"
                      << "FM: 0=none 1=Bodek-Ritchie 2=Accardi-SVG 3=Accardi-CS 4=hard-sphere 5=Wiringa\n";
            exit(0);
        } else {
            throw std::invalid_argument(std::string("Unknown argument: ") + argv[i]);
        }
    }

    if (cfg.nevent <= 0)    throw std::invalid_argument("--nevent required");
    if (cfg.output.empty()) throw std::invalid_argument("--output required");
    if (cfg.target < 0 || cfg.target > 15) throw std::invalid_argument("--target 0-15 required");
    if (cfg.nkin <= 0)      throw std::invalid_argument("--nkin required");
    if (cfg.e0 <= 0.0f)     throw std::invalid_argument("--e0 required");

    return cfg;
}

// ============================================================================
// LUND output writer (C++)
// ============================================================================

class LundWriter {
    FILE* file_;  // C FILE* kept for exact fprintf format compatibility with golden files
public:
    explicit LundWriter(const std::string& filename)
        : file_(fopen(filename.c_str(), "w")) {
        if (!file_)
            throw std::runtime_error(std::format("Cannot open output file: {}", filename));
    }

    // Non-copyable, non-movable (owns FILE*)
    LundWriter(const LundWriter&) = delete;
    LundWriter& operator=(const LundWriter&) = delete;

    void write_event(const farm::EventData& evt, int iA, int iZ, float E0,
                      int nucleon, float vz) {
        auto nb = static_cast<int>(evt.particles.size());
        fprintf(file_, "%12d%12d%12d%12d%12d%12d%12.7f%17d%12d%12.8f\n",
                nb, iA, iZ, 0, 0, 11, E0, nucleon, 1, 1.0f);
        for (int l = 0; l < nb; l++) {
            const auto& p = evt.particles[l];
            fprintf(file_, "%12d%12d%12d%12d%12d%12d%16.8E%16.8E%16.8E%16.8E%16.8E%15d%12d%12.8f\n",
                    l + 1, p.charge, 1, p.id, p.mother_id, 0,
                    p.px, p.py, p.pz, p.E, p.m, 0, 0, vz);
        }
    }

    ~LundWriter() { if (file_) fclose(file_); }
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
    farm::KinematicsResult kin_result_{};
    farm::FermiMotionState fm_state_{};
    farm::FMTableState fm_table_state_{};
    farm::AccardiFMState accardi_state_{};
    farm::DensityTable density_;
    farm::QuenchingEngine quenching_;
    std::mt19937 rng_;
    std::uniform_real_distribution<float> uniform_{0.0f, 1.0f};

public:
    explicit Simulation(const SimConfig& cfg)
        : cfg_(cfg), nkin_(cfg.nkin),
          rng_(cfg.seed >= 0 ? static_cast<unsigned>(cfg.seed) : std::random_device{}()) {}

    void initialize() {
        // Skip frequent PYTHIA re-init for targets without Fermi motion
        if (cfg_.target == 0 || cfg_.iFM == 0)
            nkin_ = cfg_.nevent + 1;

        // Nuclear parameters (C++)
        auto [Z, A, rFM] = farm::get_nuclear_params(cfg_.target);
        iZ_ = Z;
        iA_ = A;
        rFM_ = rFM;

        // Seed PYTHIA RNG
        if (cfg_.seed >= 0)
            pythia_.seed(cfg_.seed);

        // Load Fermi motion tables (C++)
        if (rFM_ != 0.0f) {
            if (cfg_.iFM == 5) {
                farm::gen_rw_table(fm_state_, iZ_, iA_, cfg_.FMlimit);
            } else if (cfg_.iFM == 2 || cfg_.iFM == 3) {
                int irho = cfg_.iFM - 1;
                farm::GenFMtable(irho, iZ_, iA_, cfg_.iFM, cfg_.target,
                                  cfg_.FMlimit, fm_table_state_, accardi_state_);
                // Copy Accardi table into FermiMotionState for sample_fermi_motion
                fm_state_.step_size = fm_table_state_.step_size;
                std::copy(fm_table_state_.FM_table.begin(),
                          fm_table_state_.FM_table.begin() + fm_table_state_.FMnb,
                          fm_state_.FM_table.begin());
            }
        }

        // Generate nuclear density table (C++)
        if (rFM_ != 0.0f)
            density_.generate(iZ_, iA_, cfg_.iDens);

        // Configure quenching engine
        quenching_ = farm::QuenchingEngine(farm::QuenchingConfig{
            .iqw = cfg_.iqw, .iqg = cfg_.iqg, .iEg = cfg_.iEg, .iPtF = cfg_.iPtF,
            .qhat = cfg_.qhat, .ehat = cfg_.ehat
        });
    }

    static constexpr int PYTHIA_REINIT_INTERVAL = 500000;
    static constexpr int PROGRESS_REPORT_INTERVAL = 10000;

    void run(LundWriter& writer) {
        int nkin_counter = nkin_;

        for (int ievent = 0; ievent < cfg_.nevent; ievent++) {
            // Re-initialize kinematics + PYTHIA when needed
            bool reinit = (ievent == 0)
                || (nkin_counter == nkin_)
                || (ievent == cfg_.nevent * iZ_ / iA_)
                || (ievent % PYTHIA_REINIT_INTERVAL == 0);

            if (reinit) {
                nkin_counter = 0;
                init_event_kinematics(ievent);
            }

            if (ievent % PROGRESS_REPORT_INTERVAL == 0)
                std::cout << std::format(" {} events processed\n", ievent);

            // Generate and process event
            auto evt = process_event();

            // Write LUND output (C++)
            float vz = farm::vertex_z(cfg_.target, uniform_(rng_), uniform_(rng_));
            writer.write_event(evt, iA_, iZ_, cfg_.e0, kin_result_.nucleon, vz);

            nkin_counter++;
        }
    }

    void print_stats() const {
        std::cout << std::format(" X sec 99 =    {:.16E}\n", pythia_.xsec(99));

        if (quenching_.QW_nb > 0)
            std::cout << std::format(" q hat =    {:E}\n", quenching_.QW_qhat / quenching_.QW_nb);
    }

private:
    void init_event_kinematics(int ievent) {
        // Fermi motion + Lorentz boost (C++)
        kin_result_ = farm::setup_kinematics(
            ievent, cfg_.nevent, iZ_, iA_,
            cfg_.iFM, rFM_, cfg_.FMlimit, cfg_.e0,
            fm_state_, rng_);

        // PYTHIA configuration + initialization (C++)
        pythia_.configure_dis();
        if (cfg_.iQuenching != 0)
            pythia_.disable_fragmentation();

        bool use_proton = (cfg_.iIso == 0)
            ? (ievent < cfg_.nevent * iZ_ / iA_)
            : (ievent < cfg_.nevent / 2);

        if (use_proton) {
            pythia_.init_proton(static_cast<double>(kin_result_.beam_energy));
        } else {
            std::cout << std::format(" X sec 99 =    {:.16E}\n", pythia_.xsec(99));
            pythia_.init_neutron(static_cast<double>(kin_result_.beam_energy));
        }
    }

    farm::EventData process_event() {
        // 1. PYTHIA event generation (C++)
        pythia_.generate_event();

        // 2. Inverse Lorentz boost back to lab frame (C++)
        if (cfg_.iFM != 0) {
            const auto& bp = kin_result_.boost;
            farm::boost_all_back(pythia_, bp.BB1, bp.B1x, bp.B1y,
                                  bp.B1z, bp.Thi, bp.Phi);
        }

        // 3. Quenching (C++)
        if (cfg_.iQuenching != 0 && cfg_.target > 1) {
            auto [ix, iy, iz] = farm::sample_interaction_position(density_, rng_);
            quenching_.apply(pythia_, ix, iy, iz,
                              density_.density_table, density_.step_size);
        }

        // 4. Fragmentation (C++ PYTHIA control)
        if (cfg_.iQuenching != 0) {
            pythia_.enable_fragmentation();
            pythia_.fragment();
            pythia_.disable_fragmentation();
        }

        // 5. Spectators (C++)
        if (cfg_.iNS == 1 && cfg_.target >= 1 && cfg_.target <= 4)
            farm::add_spectator(pythia_, cfg_.target, kin_result_.nucleon,
                                kin_result_.nuc_theta, kin_result_.nuc_phi,
                                kin_result_.nuc_momentum, kin_result_.FMintact,
                                uniform_(rng_));

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

    std::cout << " Monte Carlo Simulation Parameters:\n"
              << std::format("   Number of events:  {}\n", cfg.nevent)
              << std::format("   Lund output file:  {}\n", cfg.output)
              << std::format("   Target type:       {}\n", SimConfig::target_name(cfg.target))
              << std::format("   Nkin value:        {}\n", cfg.nkin)
              << std::format("   Electron energy:   {:.7f} GeV\n", cfg.e0)
              << std::format("   Fermi motion:      {}\n", cfg.iFM);
    if (cfg.seed >= 0)
        std::cout << std::format("   Random seed:       {}\n", cfg.seed);
    else
        std::cout << "   Random seed:       time-based\n";
    std::cout << "\n" << std::flush;

    auto t1 = std::chrono::steady_clock::now();

    // Initialize simulation
    Simulation sim(cfg);
    sim.initialize();

    // Open output (RAII — throws on failure, auto-closes on destruction)
    LundWriter writer(cfg.output);

    // Run event loop
    sim.run(writer);

    // Statistics
    sim.print_stats();
    auto t2 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t2 - t1).count();
    std::cout << std::format(" {} events in {:.2f} s\n", cfg.nevent, elapsed);

    return 0;
}

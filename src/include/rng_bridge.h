#ifndef FARM_RNG_BRIDGE_H
#define FARM_RNG_BRIDGE_H

#include <random>

// ============================================================================
// C++ RNG accessible from Fortran via extern "C"
// Replaces CERNLIB ranf() throughout the codebase.
// ============================================================================

namespace farm {

// Global RNG state (set once at initialization)
inline std::mt19937& global_rng() {
    static std::mt19937 rng;
    return rng;
}

inline void seed_global_rng(unsigned seed) {
    global_rng().seed(seed);
}

inline float random_uniform() {
    static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    return dist(global_rng());
}

} // namespace farm

// C interface callable from Fortran (replaces CERNLIB ranf)
extern "C" {
    inline float farm_random() {
        return farm::random_uniform();
    }

    inline void farm_seed_rng(int seed) {
        farm::seed_global_rng(static_cast<unsigned>(seed));
    }
}

#endif // FARM_RNG_BRIDGE_H

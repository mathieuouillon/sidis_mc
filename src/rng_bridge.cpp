#include <random>

// ============================================================================
// C++ RNG accessible from Fortran via extern "C"
// Replaces CERNLIB ranf() and ranset()
// ============================================================================

static std::mt19937 g_rng;
static std::uniform_real_distribution<float> g_dist(0.0f, 1.0f);

extern "C" {

// Called from Fortran as: result = farm_random()
// Replaces: ranf(0)
float farm_random_() {
    return g_dist(g_rng);
}

// Called from Fortran as: call farm_seed_rng(seed)
// Replaces: call ranset(seed)
void farm_seed_rng_(int* seed) {
    g_rng.seed(static_cast<unsigned>(*seed));
}

}

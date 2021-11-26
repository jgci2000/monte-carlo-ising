#include <chrono>
#include <cmath>

#include <omp.h>

#include "system.h"
#include "wl.h"
#include "rng.h"

int main(int argv, char **argc) {
    omp_set_num_threads(4);

    std::printf("Executing on %d threads \n", omp_get_num_threads());

    #pragma omp parallel
    {
        RNG rng(20 * omp_get_thread_num());
        System ising_lattice(4, 2, "SS");
        WL wang_landau(exp(1), 1 + pow(10, -8), 0.9, rng, ising_lattice);

        wang_landau.simulate(1000, omp_get_thread_num(), false);
        wang_landau.write_to_file(std::to_string(omp_get_thread_num()) + "_JDOS.txt", false);
    }
    
    return 0;
}

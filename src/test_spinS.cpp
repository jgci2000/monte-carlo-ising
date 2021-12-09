#include <chrono>
#include <cmath>

#include "system.h"
#include "wl.h"
#include "fss.h"
#include "rng.h"

int main(int argv, char **argc) {
    RNG rng(20);
    System ising_lattice(4, 2, "SS");

    double f_init = exp(1.0);
    double f_final = 1.0 + pow(10.0, -8);
    double flatness = 0.90;
    WL wang_landau(f_init, f_final, flatness, rng, ising_lattice);

    wang_landau.simulate(1000, 0, true);
    wang_landau.print_JDOS();

    // long long REP = 1000;
    // int skip = 1;
    // FSS fss(REP, skip, rng, ising_lattice);

    // fss.simulate(0, true);
    

    return 0;
}

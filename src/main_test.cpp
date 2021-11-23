
#include <iostream>
// #include <cmath>

#include "system.h"
#include "wl.h"
#include "rng.h"

int main(int argv, char **argc) {

    RNG rng(20);

    System ising_lattice(4, 2, "SS");
    ising_lattice.init_spins_random(rng);

    WL wang_landau(exp(1), 1 + pow(10, -8), 0.9, rng, ising_lattice);

    wang_landau.simulate(10000);

    return 0;
}

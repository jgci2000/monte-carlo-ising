#include <chrono>
#include <cmath>

#include "system.h"
#include "fss.h"
#include "rng.h"

int main(int argv, char **argc) {
    RNG rng(20);
    System ising_lattice(8, 2, "SS");

    FSS fss(10000, ising_lattice.N_atm, rng, ising_lattice);

    fss.simulate(0, true);
    fss.print_JDOS();
    
    return 0;
}

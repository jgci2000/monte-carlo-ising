#include <chrono>
#include <cmath>

#include "system.h"
#include "fss.h"
#include "rng.h"

int main(int argv, char **argc) {
    RNG rng(20);
    System ising_lattice(4, 2, "SS");

    FSS fss(10000, ising_lattice.N_atm, rng, ising_lattice);

    fss.simulate(0, true);
    fss.print_JDOS();
    fss.write_to_file("JDOS.txt");
    
    return 0;
}
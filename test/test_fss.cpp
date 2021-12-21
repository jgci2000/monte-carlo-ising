#include <chrono>
#include <cmath>

#include "../src/system.h"
#include "../src/fss.h"
#include "../src/rng.h"

int main(int argv, char **argc) {
    RNG rng(20);
    System ising_lattice(4, 2, "SS", "../");

    FSS fss(1000, 1, rng, ising_lattice);

    fss.simulate(0, true);
    fss.write_to_file("JDOS.txt", "./", true);
    
    return 0;
}

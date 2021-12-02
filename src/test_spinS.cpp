#include <chrono>
#include <cmath>

#include "system.h"
#include "wl.h"
#include "rng.h"

int main(int argv, char **argc) {
    RNG rng(20);
    System ising_lattice(4, 4, "SS");

    std::printf("S: %f \n", ising_lattice.S);
    std::printf("max_E: %d \n", ising_lattice.max_E);
    std::printf("max_M: %d \n", ising_lattice.max_M);

    std::printf(" spin_values: [ ");
    for (int i = 0; i < ising_lattice.Sz; ++i)
        std::printf("%d, ", ising_lattice.spins_values[i]);
    std::printf(" ] \n");

    return 0;
}

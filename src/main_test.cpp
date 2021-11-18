
#include <iostream>
// #include <cmath>

#include "system.h"
#include "wl.h"
#include "splimix64.h"

int main(int argv, char **argc) {

    // splitmix64_seed(0);

    System ising_lattice(4, 2, "SS");
    // ising_lattice.init_spins_random(&splitmix64);
    ising_lattice.init_spins_max_M();

    std::printf("original add: %p \n", &(ising_lattice));

    // std::printf("dim: %d;L: %d; Sz: %d; N_atm: %d \n", 
    // ising_lattice.dim, ising_lattice.L, ising_lattice.Sz, 
    // ising_lattice.N_atm);
    // std::printf("Lattice = %s \n", ising_lattice.lattice.data());

    // std::printf("Energies: [ ");
    // for (auto const& x : ising_lattice.energies)
    //     std::printf("%d, ", x);
    // std::printf("] \n");

    // std::printf("Magnetizations: [ ");
    // for (auto const& x : ising_lattice.magnetizations)
    //     std::printf("%d, ", x);
    // std::printf("] \n");

    // std::printf("Spins: [ ");
    // for (int i = 0; i < ising_lattice.N_atm; ++i)
    //     std::printf("%d, ", ising_lattice.spins_vector[i]);
    // std::printf(" ]\n");

    WL wang_landau(exp(1), 1 + pow(10, -8), 0.9, 10, ising_lattice);
    std::printf("done2!\n");

    // std::printf("f_vals: [ ");
    // for (int i = 0; i < wang_landau.n_f_vals; ++i) 
    //     std::printf("%f, ", wang_landau.f_vals[i]);
    // std::printf(" ]\n");

    return 0;
}


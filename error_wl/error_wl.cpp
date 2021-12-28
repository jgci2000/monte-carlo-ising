#include <iostream>
#include <cmath>

#include <omp.h>

#include "../src/wl.h"
#include "../src/rng.h"
#include "../src/system.h"

int main(int argc, char **argv) {
    int run_max = 1000;
    int L = 4;
    int f_final_exp = 10;
    int flatness_vals[] = {70, 80, 90, 95, 99}; int n = 5;
    double f_init = exp(1);

    #pragma omp parallel for 
    for (int f_exp = 1; f_exp < f_final_exp; f_exp++) {
        double f_final = 1.0 + pow(10, -f_exp);

        for (int flatness_idx = 0; flatness_idx < n; flatness_idx++) {
            double flatness = flatness_vals[flatness_idx] / 100.0;

            for (int run = 0; run < run_max; run++) {
                std::string path = "data/f_" + std::to_string(f_exp) + "/flatness_" + std::to_string(flatness_vals[flatness_idx]) + "/";
                std::string name = std::to_string(run) + "_JDOS.txt";

                RNG rng((unsigned) time(NULL));
                System ising_lattice(L, 2, "SS", "../");

                WL wang_landau(f_init, f_final, flatness, rng, ising_lattice);
                wang_landau.simulate(1000, run, false);
                wang_landau.write_to_file(name, path, true);
            }
        }        
    }

    return 0;
}

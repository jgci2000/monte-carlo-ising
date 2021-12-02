#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include <omp.h>

#include "../src/wl.h"
#include "../src/rng.h"
#include "../src/system.h"

int main(int argc, char **argv) {

    int L, f_init_exp, f_final_exp, run;
    double flatness;

    run = atoi(argv[1]);
    L = atoi(argv[2]);
    f_init_exp = atoi(argv[3]);
    f_final_exp = atoi(argv[4]);
    flatness = atoi(argv[5]) / 100.0;

    std::vector<double> error;
    std::vector<double> f_vals;
    
    double f_init = exp(1);
    #pragma omp parallel for 
    for (int f_exp = f_init_exp; f_exp < f_final_exp; f_exp++) {
        int status = system("mkdir data");
        std::string dir = "data/" + std::to_string(atoi(argv[5])) + "/" + std::to_string(f_exp);
        double f_final = 1.0 + pow(10, -f_exp);

        for (int i = 0; i < run; ++i) {
            std::string file_name = std::to_string(i) + "_JDOS.txt";

            RNG rng((unsigned) time(NULL));
            System ising_lattice(L, 2, "SS", "../");

            WL wang_landau(f_init, f_final, flatness, rng, ising_lattice);
            wang_landau.simulate(1000, i, false);
            wang_landau.write_to_file(file_name, dir, false);
        }
    }




    return 0;
}

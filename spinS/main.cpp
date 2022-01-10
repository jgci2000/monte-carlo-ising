#include <iostream>
#include <omp.h>

#include "../src/fss.h"
#include "../src/rng.h"
#include "../src/system.h"

int main(int argc, char **argv) {
    int L = 4;
    std::string lattice = "SS";
    int Sz_vals[] = {2, 3, 4, 5};    
    int n_Sz = 4;

    int REP = 10000; int REP_exp = 4;
    int skip = L * L / 4;
    
    #pragma omp parallel for
    for (int i = 0; i < n_Sz; i++) {
        double Sz = Sz_vals[i];
        std::string filename = "JDOS_L" + std::to_string(L) + "_" + lattice +  "_Sz_" + std::to_string(Sz) + "_R1E" + std::to_string(REP_exp) + ".txt";
        std::string path = "data/";

        RNG rng((unsigned) time(NULL));
        System ising(L, Sz, lattice, "../");

        FSS fss(REP, skip, rng, ising);;
        fss.simulate();
        fss.write_to_file(filename, path);
    }

    return 0;
}



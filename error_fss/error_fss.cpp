#include <iostream>
#include <cmath>

#include <omp.h>

#include "../src/fss.h"
#include "../src/rng.h"
#include "../src/system.h"

int main(int argc, char** argv) {
    int run_max = 1000;
    int L = 4; int N_atm = L * L;
    int exp_REP_final = 7;
    int skip_vals[] = {1, N_atm/4, N_atm/2, N_atm, 2*N_atm, 4*N_atm}; int n = 6;
    std::string skip_str[] = {"0", "0.25", "0.5", "1", "2", "4"};
    
    omp_set_num_threads(6);
    #pragma omp parallel for
    for (int REP_exp = 2; REP_exp <= exp_REP_final; REP_exp++) {
        long long REP = pow(10, REP_exp);
        
        for (int skip_idx = 0; skip_idx < n; skip_idx++) {
            int skip = skip_vals[skip_idx];

            for (int run = 0; run < run_max; run++) {
                std::string path = "data/REP_" + std::to_string(REP_exp) + "/skip_" + skip_str[skip_idx] + "/";
                std::string name = std::to_string(run) + "_JDOS.txt";

                RNG rng((unsigned) time(NULL));
                System ising(L, 2, "SS", "../");

                FSS fss(REP, skip, rng, ising);
                fss.simulate(run, false);
                fss.write_to_file(name, path, true);
            }
        }
    }    


    return 0;
}

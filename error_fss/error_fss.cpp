#include <iostream>
#include <cmath>

#include <omp.h>

#include "../src/fss.h"
#include "../src/rng.h"
#include "../src/system.h"

int main(int argc, char** argv) {
    int L, run_max;
    int exp_REP_final;

    run_max = atoi(argv[1]);
    L = atoi(argv[2]);
    exp_REP_final = atoi(argv[3]);
    
    #pragma omp parallel for
    for (int REP_exp = 2; REP_exp <= exp_REP_final; REP_exp++) {
        long long REP = pow(10, REP_exp);
        
        for (int skip_exp = -2; skip_exp <= 2; skip_exp++) {
            int skip = L * L * pow(2, skip_exp);

            for (int run = 1; run <= run_max; run++) {
                std::string path = "data/REP_" + std::to_string(REP_exp) + "/skip_" + std::to_string(skip_exp) + "/";
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

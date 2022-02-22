#include <iostream>
#include <cmath>
#include <fstream>

#include "../src/fss.h"
#include "../src/rng.h"
#include "../src/system.h"

int main(int argc, char** argv) {
    int L = atoi(argv[1]); 
    int N_atm = L * L;
    long REP = 10000;
    int skip = N_atm / 4;

    std::string path = "data/";
    std::string name = "JDOS_L" + std::to_string(L) + "_SS_Sz2_R1E4_skip" + std::to_string(skip) + ".txt";

    RNG rng((unsigned) time(NULL));
    System ising(L, 2, "SS", "../");

    FSS fss(REP, skip, rng, ising);
    fss.simulate(0, true);
    fss.write_to_file(name, path, true);

    std::ofstream run_time_file("run_times.txt", std::ios::app);
    if (run_time_file.is_open()) {
        run_time_file << "L" + std::to_string(L) + "_SS_Sz2_R1E4_skip" + std::to_string(skip) + ": " + std::to_string(fss.run_time) + "\n"; 
    }
    run_time_file.close();
    
    return 0;
}

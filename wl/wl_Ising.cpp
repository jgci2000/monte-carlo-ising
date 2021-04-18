//
// Wang Landau sampling for the Ising 1/2 Model 
// João Inácio, Apr. 18th, 2021
//
// This version is single core
//

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdint.h>
#include <chrono>
#include <cmath>
#include <array>
#include <vector>
#include <map>
#include <filesystem>

#include "WL_Functions.h"
#include "splimix64.h"
#include "xoshiro256++.h"

using std::cout;
using std::endl;
using std::vector;
using std::array;
using std::string;
using std::map;
using std::to_string;

#define ll          long long
#define ld          long double


// Seed for the RNG. If SEED equal 0, use random seed.
#define SEED        0

// Size of the Ising Lattice
#define L_LATTICE   4
// LATTICE_NUM -> 1 - SS; 2 - SC; 3 - BCC; 4 - FCC; 5 - HCP; 6 - Hex 
#define LATTICE_NUM 1

// Output location
#define SAVE_DIR(lattice, L, f_final)    "./data/" + lattice + "/L" + to_string(L) + "/" + to_string(f_final) + "/"


int main(int argc, char **argv)
{
    // Set the seed for xoshiro256++

    uint64_t seed = SEED;
    if (seed == 0)
    {
        srand((unsigned) time(NULL));
        seed = rand();
    }

    splitmix64_seed(seed);
    for (int i = 0; i < 4; i++)
        s[i] = splitmix64();
    
    // Initialize Ising and set parameters for FSS computations
    
    system_info system = get_system(L_LATTICE, LATTICE_NUM);

    int L = system.L;
    string lattice = system.lattice;
    int dim = system.dim;
    int NN = system.NN;
    int N_atm = system.N_atm;

    int max_E = (1.0 / 2.0) * NN * N_atm;
    int max_M = N_atm;

    int NE = 1 + (max_E / 2);
    int NM = N_atm + 1;

    map<int, int> energies = create_map(- max_E, max_E, 4);
    map<int, int> magnetizations = create_map(- max_M, max_M, 2);

    int run;
    double f = exp(1);
    double f_final;
    double flatness;

    switch (argc)
    {
        case 2:
            f_final = 1 + pow(10, - atoi(argv[1]));
            run = 0;
            flatness = 0.9;
            break;
        case 3:
            f_final = 1 + pow(10, - atoi(argv[1]));
            run = atoi(argv[2]);
            flatness = 0.9;
            break;
        case 4:
            f_final = 1 + pow(10, - atoi(argv[1]));
            flatness = atof(argv[2]) / 100;
            run = atoi(argv[3]);
            break;
        default:
            cout << "No parameters selected, reseting do default." << endl;
            run = 0;
            flatness = 0.9;
            f_final = f_final = 1 + pow(10, - 8);
            break;
    }

    string NN_table_file_name = "./neighbour_tables/neighbour_table_" + to_string(dim) + "D_" + lattice + 
    "_" + to_string(NN) + "NN_L" + to_string(L) + ".txt";
    string norm_factor_file_name = "./coefficients/coefficients_" + to_string(N_atm) + "d2.txt";
    string save_file = to_string(run) + "_JDOS_WL_Ising_" + to_string(dim) + "D_" + lattice + "_L" + to_string(L) + "_f" + 
    to_string((int) - log10(f_final - 1)) + "_flatness" + to_string((int) (flatness * 100));

    // Initialize vectors and read files

    int *spins_vector = new int[N_atm];
    int *NN_table = new int[N_atm * NN];
    ld *norm_factor = new ld[NM];
    
    read_NN_talbe(NN_table_file_name, NN_table);
    read_norm_factor(norm_factor_file_name, norm_factor);

    ld *ln_JDOS = new ld[NE * NM];
    ld *JDOS = new ld[NE * NM];
    ll *hist = new ll[NE * NM];

    for (int i = 0; i < NE * NM; i++)
    {
        JDOS[i] = 0;
        ln_JDOS[i] = 0;
        hist[i] = 0;
    }

    ll mc_sweep = 0;

    // Random spins config

    for (int i = 0; i < N_atm; i++) 
    {
        if ((rand_xoshiro256pp() % 2) + 1 == 1) 
            spins_vector[i] = + 1;
        else spins_vector[i] = - 1;
    }

    int E_config = 0;
    int M_config = 0;
    for (int i = 0; i < N_atm; i++)
    {
        for (int a = 0; a < NN; a++)
            E_config += - spins_vector[i] * spins_vector[NN_table[i * NN + a]];

        M_config += spins_vector[i];
    }
    E_config = E_config / 2;

    int idx_E_config = energies[E_config];
    int idx_M_config = magnetizations[M_config];

    // Start measuring time

    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;

    now = time(0);
    t = ctime(&now); t.pop_back();

    string console_output = "run: " + to_string(run) + " | L: " + to_string(L) + " | f_final: 1+1E" + to_string((int) log10(f_final - 1)) + 
    " | flatness: " + to_string((int) (flatness * 100)) + " | dim: " + to_string(dim) + "D | lattie: " + lattice;
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();

    // Wang Landau Sampling
    // Main loop

    std::chrono::_V2::steady_clock::time_point loop_start;
    
    while(f > f_final)
    {
        if (mc_sweep == 0)
            loop_start = std::chrono::steady_clock::now();
        
        for (int idx = 0; idx < N_atm; idx++) 
        {
            int flip_idx = (rand_xoshiro256pp() % N_atm);

            int delta_E = 0;
            for (int a = 0; a < NN; a++)
                delta_E += - spins_vector[flip_idx] * spins_vector[NN_table[flip_idx * NN + a]];
            
            int new_E_config = E_config - 2 * delta_E;
            int new_M_config  = M_config - 2 * spins_vector[flip_idx];
            int new_idx_E_config = energies[new_E_config];
            int new_idx_M_config = magnetizations[new_M_config];
            
            ld ratio = exp(ln_JDOS[idx_E_config * NM + idx_M_config] - ln_JDOS[new_idx_E_config * NM + new_idx_M_config]);

            if (ratio >= 1 || ((ld) rand_xoshiro256pp() / (ld) UINT64_MAX) < ratio)
            {
                spins_vector[flip_idx] = - spins_vector[flip_idx];
                
                E_config = new_E_config;
                idx_E_config = new_idx_E_config;
                M_config = new_M_config;
                idx_M_config = new_idx_M_config;
            }

            hist[idx_E_config * NM + idx_M_config]++;
            ln_JDOS[idx_E_config * NM + idx_M_config] += log(f);
        }

        mc_sweep++;

        if (mc_sweep % 10000 == 0) 
        {
            long double avg_h = average_hist(hist, NE * NM);
            int min_h = min_hist(hist, NE * NM);
            
            if (min_h >= avg_h * flatness)
            {
                auto loop_end = std::chrono::steady_clock::now();
                double loop_dur = (double) (std::chrono::duration_cast<std::chrono::microseconds> (loop_end - loop_start).count()) * pow(10, -6);

                now = time(0);
                t = ctime(&now); t.pop_back();
                
                string console_output = t + " | f: 1+1E" + to_string(log10(f - 1)) + "/1+1E" + to_string((int) log10(f_final - 1)) + 
                " | sweeps: " + to_string(mc_sweep) + " | flat time: " + to_string(loop_dur) + "s";
                string data_line = "1+1E" + to_string(log10(f - 1)) + " 1+1E" + to_string((int) log10(f_final - 1)) + " " + 
                to_string(loop_dur) + " " + to_string(mc_sweep) + " " + to_string(min_h) + " " + to_string(avg_h);

                console_log.push_back(console_output);
                data.push_back(data_line);

                cout << console_output << endl;

                f = sqrt(f);
                mc_sweep = 0;

                for (int i = 0; i < NE * NM; i++)
                    hist[i] = 0;
            }
        }
    }

    // Normalize JDOS

    for (int q = 0; q < NM; q++)
    {
        int first_idx;
        for (int i = 0; i < NE; i++)
        {
            if (ln_JDOS[i * NE + q] > 0)
            {
                first_idx = i;
                break;
            }
        }

        ld temp = 0;
        for (int i = 0; i < NE; i++)
            if (ln_JDOS[i * NM + q] > 0)
                temp += exp(ln_JDOS[i * NM + q] - ln_JDOS[first_idx * NM + q]);

        ld sum_ln_JDOS = ln_JDOS[first_idx * NM + q] + log(temp);
        
        for (int i = 0; i < NE; i++)
            if (ln_JDOS[i * NM + q] > 0)
                JDOS[i * NM + q] = exp(ln_JDOS[i * NM + q] + log(norm_factor[q]) - sum_ln_JDOS);  
    }

    // Stop mesuring time

    auto method_end = std::chrono::steady_clock::now();
    double runtime = (double) (std::chrono::duration_cast<std::chrono::microseconds> (method_end - method_start).count()) * pow(10, -6);
    now = time(0);
    t = ctime(&now); t.pop_back();

    cout << endl;
    cout << "Runtime: " << std::setw(8) << runtime << " seconds." << endl;
    cout << "Simulation ended at: " << t << endl;

    // Write JDOS to file

    namespace fs = std::filesystem;
    fs::create_directories(SAVE_DIR(lattice, L, (int) - log10(f_final - 1)));

    std::ofstream file1((string) SAVE_DIR(lattice, L, (int) - log10(f_final - 1)) + save_file + ".txt");
    for (int i = 0; i < NE; i++) 
    {
        for (int j = 0; j < NM; j++) 
            file1 << JDOS[i * NM + j] << " ";
        file1 << "\n";
    }
    file1.close();

    std::ofstream file2((string) SAVE_DIR(lattice, L, (int) - log10(f_final - 1)) + save_file + "_data.txt");
    file2 << "f f_max loop_dur mc_sweeps min_h avg_h \n"; 
    for (int i = 0; i < data.size(); i++)
        file2 << data[i] << "\n";
    file2 << runtime << "\n";
    file2.close();

    std::ofstream file3((string) SAVE_DIR(lattice, L, (int) - log10(f_final - 1)) + save_file + "_console_logs.txt");
    for (int i = 0; i < console_log.size(); i++)
        file3 << console_log.at(i) << "\n";
    file3.close();

    delete[] JDOS, ln_JDOS, hist;
    delete[] spins_vector;
    delete[] NN_table, norm_factor;

    return 0;
}

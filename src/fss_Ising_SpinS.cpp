//
// Flat Scan Sampling for the Ising SpinS Model 
// João Inácio, May 8th, 2021
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

#include "../include/functions.h"
#include "../include/splimix64.h"
#include "../include/xoshiro256++.h"

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

// Quantum number of the spin of the system
// This is used to compute the SZ of the particles
// SZ = 2*S + 1
#define S           (double) 1
// Size of the Ising Lattice
#define L_LATTICE   4
// LATTICE_NUM -> 1 - SS; 2 - SC; 3 - BCC; 4 - FCC; 5 - HCP; 6 - Hex 
#define LATTICE_NUM 1

// Output location
#define SAVE_DIR(lattice, L, log_REP)    "./data/fss_spinS/" + lattice + "/L" + to_string(L) + "/" + to_string(log_REP) + "/"


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
    
    system_info system = get_system(L_LATTICE, LATTICE_NUM, S);

    int L = system.L;
    string lattice = system.lattice;
    int dim = system.dim;
    int NN = system.NN;
    int N_atm = system.N_atm;
    int SZ = system.SZ;

    const int max_E = 4 * S * S * N_atm * NN / 2;
    const int max_M = 2 * S * N_atm;

    int NE = 1 + (max_E / 2);
    int NM = max_M + 1;

    map<int, int> energies = create_map(- max_E, max_E, 4);
    map<int, int> magnetizations = create_map(- max_M, max_M, 2);
    
    vector<int> spinZ;
    for (int i = 0; i < SZ; i++)
        spinZ.push_back(- 2 * S + 2 * i);
    
    int q_max = (NM + 1) / 2 - 2;
    if (NM % 2 == 0)
        q_max = NM / 2 - 3;

    int run;
    int skip;
    ll REP;

    switch (argc)
    {
        case 2:
            REP = pow(10, atoi(argv[1]));
            run = 0;
            skip = N_atm;
            break;
        case 3:
            REP = pow(10, atoi(argv[1]));
            run = atoi(argv[2]);
            skip = N_atm;
            break;
        case 4:
            REP = pow(10, atoi(argv[1]));
            skip = atoi(argv[2]);
            run = atoi(argv[3]);
            break;
        default:
            cout << "No parameters selected, reseting do default." << endl;
            run = 0;
            skip = N_atm;
            REP = pow(10, 4);
            break;
    }

    string NN_table_file_name = "./neighbour_tables/neighbour_table_" + to_string(dim) + "D_" + 
    lattice + "_" + to_string(NN) + "NN_L" + to_string(L) + ".txt";
    string norm_factor_file_name = "./coefficients/coefficients_" + to_string(N_atm) + "d" + to_string(SZ) + ".txt";
    string Npos_file_name = "./sum_npos/sum_configs_Npos" + std::to_string(SZ) + "_N_atm" + std::to_string(N_atm) + ".txt";
    string save_file = to_string(run) + "_JDOS_FSS_Ising_SpinS_" + to_string(dim) + "D_" + lattice + "_SZ" + to_string(SZ) + 
    "_L" + to_string(L) + "_REP_1E" + to_string((int) log10(REP)) + "_skip_" + to_string(skip);

    // Initialize vectors and read files

    int *spins_vector = new int[N_atm];
    int *new_spins_vector = new int[N_atm];
    int *NN_table = new int[N_atm * NN];
    ld *norm_factor = new ld[NM];
    ll *hist;
    ll *hist_E_selected;

    int **Npos;
    vector<int> line_size_Npos;
    
    read_NN_talbe(NN_table_file_name, NN_table);
    read_norm_factor(norm_factor_file_name, norm_factor);

    string line;
    std::ifstream sum_npos_file(Npos_file_name);
    int size = 0;

    if (sum_npos_file.is_open())
    {
        int i = 0;

        std::getline(sum_npos_file, line);
        line_size_Npos = split_int(line, ' ');

        Npos = new int *[line_size_Npos.size()];

        while (std::getline(sum_npos_file, line))
        {
            Npos[i] = new int[line_size_Npos.at(i) * SZ];
            vector<int> a = split_int(line, ' ');

            for (int k = 0; k < line_size_Npos.at(i) * SZ; k++)
                Npos[i][k] = a.at(k);
            i++;
        }
        sum_npos_file.close();
    }
    else
    {
        cout << "Unable to sum Nppos file. Invalid lattice size, lattice type or total spin." << endl;
    }

    ld *JDOS = new ld[NE * NM];
    ld **JDOS_M_spin = new ld *[NM];

    for (int i = 0; i < NE * NM; i++)
        JDOS[i] = 0;
    JDOS[0] = 1;
    
    for (int i = 0; i < NM; i++)
        JDOS_M_spin[i] = new ld[line_size_Npos.at(i) * NE];
    
    for (int i = 0; i < NM; i++)
        for (int j = 0; j < NE * line_size_Npos.at(i); j++)
            JDOS_M_spin[i][j] = 0;
    JDOS_M_spin[0][0] = 1;
    
    // Start measuring time

    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;

    now = time(0);
    t = ctime(&now); t.pop_back();

    string console_output = "run: " + to_string(run) + " | L: " + to_string(L) + " | REP: " + to_string(REP) + " | skip: " + 
    to_string(skip) + " | dim: " + to_string(dim) + "D | lattice: " + lattice + " | SZ: " + to_string(SZ);
    console_log.push_back(console_output);

    cout << endl;
    cout << console_output << endl;
    cout << "Starting time: " << t << endl << endl;

    auto method_start = std::chrono::steady_clock::now();

    // Flat Scan Sampling
    // Scan and Compute JDOS at q = 1

    vector<int> SPM; 
    for (int i = 0; i < N_atm; i++)
        SPM.push_back(0);

    int prev_idx_E_config = 0;
    int prev_idx_Npos = 0;

    for (int x = 1; x < SZ; x++)
    {
        for (int flip_idx = 0; flip_idx < N_atm; flip_idx++)
        {
            vector<int> SPM_tmp = SPM;

            int E_tmp1 = 0;
            for (int a = 0; a < NN; a++)
                E_tmp1 += - spinZ[SPM[flip_idx]] * spinZ[SPM[flip_idx]];
            
            int E_tmp2 = 0;
            for (int a = 0; a < NN; a++)
                E_tmp2 += - spinZ[SPM[flip_idx] + x] * spinZ[SPM[flip_idx]];
            
            int E_tmp3 = - max_E - E_tmp1 + E_tmp2;
            int idx_E_tmp3 = energies[E_tmp3];

            SPM_tmp[flip_idx] = SPM[flip_idx] + x;

            vector<int> counter; counter.assign(SZ, 0);
            for (int i = 0; i < N_atm; i++)
                for (int j = 0; j < SZ; j++)
                    if (SPM_tmp[i] == j)
                        counter[j]++;
            
            int counter2 = 0;
            int idx_Npos;
            for (idx_Npos = 0; idx_Npos < line_size_Npos.at(x) * SZ; idx_Npos += SZ)
            {
                for (int i = 0; i < SZ; i++)
                    if (counter[i] == Npos[x][i + idx_Npos])
                        counter2++;

                if (counter2 == SZ)
                    break;
                else
                    counter2 = 0;
            }
            idx_Npos /= SZ;

            JDOS_M_spin[x][idx_Npos * NE + idx_E_tmp3] += JDOS_M_spin[0][0];
        }
    }

    vector<ld> sum_JDOS_M_spin; sum_JDOS_M_spin.assign(NE, 0);
    ld sum_sum_JDOS_M_spin = 0;
    for (int i = 0; i < NE; i++)
    {
        for (int j = 0; j < line_size_Npos.at(1); j++)
            sum_JDOS_M_spin[i] += JDOS_M_spin[1][j * NE + i];
        sum_sum_JDOS_M_spin += sum_JDOS_M_spin[i];
    }

    for (int i = 0; i < NE; i++)
        JDOS[i * NM + 1] = sum_JDOS_M_spin[i] * norm_factor[1] / sum_sum_JDOS_M_spin;

    console_output = t + " | q: " + to_string(0) + "/" + to_string(q_max);
    console_log.push_back(console_output);

    cout << console_output << endl;

    // Main Loop

    for (int q = 1; q <= q_max; q++)
    {
        auto q_start = std::chrono::steady_clock::now();

        hist = new ll[NE * line_size_Npos.at(q)];
        hist_E_selected = new ll[NE * line_size_Npos.at(q)];

        for (int i = 0; i < NE * line_size_Npos.at(q); i++)
        {
            hist[i] = 0; 
            hist_E_selected[i] = 0;
        }

        // Random config at q

        for (int i = 0; i < N_atm; i++)
            spins_vector[i] = spinZ[0];
        int E_config =  - max_E;
        int idx_E_config = energies[E_config];

        for (int i = 0; i < N_atm; i++)
            SPM[i] = 0;

        for (int idx = 1; idx <= q; idx++)
        {
            vector<int> flip_list;
            for (int i = 0; i < N_atm; i++)
                if (SPM[i] + 1 < SZ)
                    flip_list.push_back(i);
            int flipped_idx = flip_list.at(rand_xoshiro256pp() % flip_list.size());

            int E_old = 0;
            for (int a = 0; a < NN; a++)
                E_old += - spins_vector[flipped_idx] * spins_vector[NN_table[flipped_idx * NN + a]];

            SPM[flipped_idx]++;
            spins_vector[flipped_idx] = spinZ[SPM[flipped_idx]];

            int E_new = 0;
            for (int a = 0; a < NN; a++)
                E_new += - spins_vector[flipped_idx] * spins_vector[NN_table[flipped_idx * NN + a]];
            
            E_config = E_config - E_old + E_new;
        }
        
        vector<int> counter; counter.assign(SZ, 0);
        for (int i = 0; i < N_atm; i++)
            for (int j = 0; j < SZ; j++)
                if (SPM[i] == j)
                    counter[j]++;
        
        int counter2 = 0;
        int idx_Npos;
        for (idx_Npos = 0; idx_Npos < line_size_Npos.at(q) * SZ; idx_Npos += SZ)
        {
            for (int i = 0; i < SZ; i++)
                if (counter[i] == Npos[q][i + idx_Npos])
                    counter2++;

            if (counter2 == SZ)
                break;
            else
                counter2 = 0;
        }
        
        idx_E_config = energies[E_config];
        int idx_Npos_config = idx_Npos / SZ;

        // Update Histograms

        hist[idx_E_config * line_size_Npos.at(q) + idx_Npos_config]++;
        hist_E_selected[idx_E_config * line_size_Npos.at(q) + idx_Npos_config]++;
        
        // Scan the first config

        for (int x = 1; x < SZ; x++)
        {
            vector<int> flip_list;
            for (int i = 0; i < N_atm; i++)
                if (SPM[i] + 1 <= SZ - x)
                    flip_list.push_back(i);

            for (int flip_idx = 0; flip_idx < flip_list.size(); flip_idx++)
            {
                vector<int> SPM_tmp = SPM;

                int E_tmp1 = 0;
                for (int a = 0; a < NN; a++)
                    E_tmp1 += - spins_vector[flip_list.at(flip_idx)] * spins_vector[NN_table[flip_list.at(flip_idx) * NN + a]];

                int E_tmp2 = 0;
                for (int a = 0; a < NN; a++)
                    E_tmp2 += - spinZ[SPM[flip_list.at(flip_idx)] + x] * spins_vector[NN_table[flip_list.at(flip_idx) * NN + a]];

                int E_tmp3 = E_config - E_tmp1 + E_tmp2;
                int idx_E_tmp3 = energies[E_tmp3];

                SPM_tmp[flip_list.at(flip_idx)] = SPM[flip_list.at(flip_idx)] + x;

                vector<int> counter; counter.assign(SZ, 0);
                for (int i = 0; i < N_atm; i++)
                    for (int j = 0; j < SZ; j++)
                        if (SPM_tmp[i] == j)
                            counter[j]++;
                
                int counter2 = 0;
                int idx_Npos;
                for (idx_Npos = 0; idx_Npos < line_size_Npos.at(q + x) * SZ; idx_Npos += SZ)
                {
                    for (int i = 0; i < SZ; i++)
                        if (counter[i] == Npos[q + x][i + idx_Npos])
                            counter2++;

                    if (counter2 == SZ)
                        break;
                    else
                        counter2 = 0;
                }
                idx_Npos /= SZ;

                JDOS_M_spin[q + x][idx_Npos * NE + idx_E_tmp3] += JDOS_M_spin[q][idx_Npos_config * NE + idx_E_config] / REP;
            }
        }

        ll k = 1;
        bool accepted = false;
        vector<int> new_SPM;

        // Where the magic happens

        while (min_hist(hist_E_selected, NE * line_size_Npos.at(q)) < REP)
        {
            // Get a new random condig at magnetization q
            
            if (!accepted)
            {
                new_SPM = SPM;
                for (int i = 0; i < N_atm; i++)
                    new_spins_vector[i] =  spins_vector[i];
            }
            int new_E_config = 0;
            int new_idx_E_config = 0;
            int new_idx_Npos_config = 0;

            // Choose a spin to flip
            
            int flipped_idx1 = rand_xoshiro256pp() % N_atm;

            int E_old = 0;
            for (int a = 0; a < NN; a++)
                E_old += - new_spins_vector[flipped_idx1] * new_spins_vector[NN_table[flipped_idx1 * NN + a]];

            vector<int> SPM_end_list;
            for (int i = 0; i < SZ; i++)
                if (i != SPM[flipped_idx1])
                    SPM_end_list.push_back(i);

            int SPM_start = SPM[flipped_idx1];
            int SPM_end = SPM_end_list.at(rand_xoshiro256pp() % SPM_end_list.size());

            new_spins_vector[flipped_idx1] = spinZ[SPM_end];
            new_SPM[flipped_idx1] = SPM_end;

            int E_new = 0;
            for (int a = 0; a < NN; a++)
                E_new += - new_spins_vector[flipped_idx1] * new_spins_vector[NN_table[flipped_idx1 * NN + a]];
            new_E_config = E_config - E_old + E_new;

            int SPM_dif = SPM_end - SPM_start;
            
            // Choose another one to flip back
            
            vector<int> flip_list;
            for (int i = 0; i < N_atm; i++)
                if (new_SPM[i] + 1 - SPM_dif >= 1 && new_SPM[i] + 1 - SPM_dif <= SZ)
                    flip_list.push_back(i);

            int flipped_idx2 = flip_list.at(rand_xoshiro256pp() % flip_list.size());

            E_old = 0;
            for (int a = 0; a < NN; a++)
                E_old += - new_spins_vector[flipped_idx2] * new_spins_vector[NN_table[flipped_idx2 * NN + a]];

            new_SPM[flipped_idx2] = new_SPM[flipped_idx2] - SPM_dif;
            new_spins_vector[flipped_idx2] = spinZ[new_SPM[flipped_idx2]];

            E_new = 0;
            for (int a = 0; a < NN; a++)
                E_new += - new_spins_vector[flipped_idx2] * new_spins_vector[NN_table[flipped_idx2 * NN + a]];
            new_E_config = new_E_config - E_old + E_new;

            vector<int> counter; counter.assign(SZ, 0);
            for (int i = 0; i < N_atm; i++)
                for (int j = 0; j < SZ; j++)
                    if (new_SPM[i] == j)
                        counter[j]++;
            
            int counter2 = 0;
            int idx_Npos;
            for (idx_Npos = 0; idx_Npos < line_size_Npos.at(q) * SZ; idx_Npos += SZ)
            {
                for (int i = 0; i < SZ; i++)
                    if (counter[i] == Npos[q][i + idx_Npos])
                        counter2++;

                if (counter2 == SZ)
                    break;
                else
                    counter2 = 0;
            }
            
            new_idx_Npos_config = idx_Npos / SZ;
            new_idx_E_config = energies[new_E_config];
            
            // Wang Landau criteria

            ld ratio = JDOS_M_spin[q][idx_Npos_config * NE + idx_E_config] / JDOS_M_spin[q][new_idx_Npos_config * NE + new_idx_E_config];

            if (ratio >= 1 || ((ld) rand_xoshiro256pp() / (ld) UINT64_MAX) < ratio || hist_E_selected[new_idx_E_config * line_size_Npos.at(q) + new_idx_Npos_config] == 0)
            {
                for (int i = 0; i < N_atm; i++)
                    spins_vector[i] = new_spins_vector[i];

                E_config = new_E_config;
                idx_E_config = new_idx_E_config;
                idx_Npos_config = new_idx_Npos_config;
                SPM = new_SPM;

                accepted = true;
            }
            else
            {
                accepted = false;
            }

            hist[idx_E_config * line_size_Npos.at(q) + idx_Npos_config]++;

            // Scan configuration

            if (hist_E_selected[idx_E_config * line_size_Npos.at(q) + idx_Npos_config] < REP && k % skip == 0 || hist_E_selected[new_idx_E_config * line_size_Npos.at(q) + new_idx_Npos_config] == 0)
            {
                for (int x = 1; x < SZ; x++)
                {
                    vector<int> flip_list;
                    for (int i = 0; i < N_atm; i++)
                        if (SPM[i] + 1 <= SZ - x)
                            flip_list.push_back(i);

                    for (int flip_idx = 0; flip_idx < flip_list.size(); flip_idx++)
                    {
                        vector<int> SPM_tmp = SPM;

                        int E_tmp1 = 0;
                        for (int a = 0; a < NN; a++)
                            E_tmp1 += - spins_vector[flip_list.at(flip_idx)] * spins_vector[NN_table[flip_list.at(flip_idx) * NN + a]];

                        int E_tmp2 = 0;
                        for (int a = 0; a < NN; a++)
                            E_tmp2 += - spinZ[SPM[flip_list.at(flip_idx)] + x] * spins_vector[NN_table[flip_list.at(flip_idx) * NN + a]];

                        int E_tmp3 = E_config - E_tmp1 + E_tmp2;
                        int idx_E_tmp3 = energies[E_tmp3];

                        SPM_tmp[flip_list.at(flip_idx)] = SPM[flip_list.at(flip_idx)] + x;

                        vector<int> counter; counter.assign(SZ, 0);
                        for (int i = 0; i < N_atm; i++)
                            for (int j = 0; j < SZ; j++)
                                if (SPM_tmp[i] == j)
                                    counter[j]++;
                        
                        int counter2 = 0;
                        int idx_Npos;
                        for (idx_Npos = 0; idx_Npos < line_size_Npos.at(q + x) * SZ; idx_Npos += SZ)
                        {
                            for (int i = 0; i < SZ; i++)
                                if (counter[i] == Npos[q + x][i + idx_Npos])
                                    counter2++;

                            if (counter2 == SZ)
                                break;
                            else
                                counter2 = 0;
                        }
                        idx_Npos /= SZ;

                        JDOS_M_spin[q + x][idx_Npos * NE + idx_E_tmp3] += JDOS_M_spin[q][idx_Npos_config * NE + idx_E_config] / REP;
                    }
                }

                hist_E_selected[idx_E_config * line_size_Npos.at(q) + idx_Npos_config]++;
            }

            k++;
        }

        // Normalize JDOS and output to console

        vector<ld> sum_JDOS_M_spin; sum_JDOS_M_spin.assign(NE, 0);
        ld sum_sum_JDOS_M_spin = 0;
        for (int i = 0; i < NE; i++)
        {
            for (int j = 0; j < line_size_Npos.at(q + 1); j++)
                sum_JDOS_M_spin[i] += JDOS_M_spin[q + 1][j * NE + i];
            sum_sum_JDOS_M_spin += sum_JDOS_M_spin[i];
        }

        for (int i = 0; i < NE * line_size_Npos.at(q + 1); i++)
            JDOS_M_spin[q + 1][i] = JDOS_M_spin[q + 1][i] * norm_factor[q + 1] / sum_sum_JDOS_M_spin;
        
        for (int i = 0; i < NE; i++)
            sum_JDOS_M_spin[i] = 0;
        sum_sum_JDOS_M_spin = 0;
        for (int i = 0; i < NE; i++)
        {
            for (int j = 0; j < line_size_Npos.at(q + 1); j++)
                sum_JDOS_M_spin[i] += JDOS_M_spin[q + 1][j * NE + i];
            sum_sum_JDOS_M_spin += sum_JDOS_M_spin[i];
        }

        for (int i = 0; i < NE; i++)
            JDOS[i * NM + q + 1] = sum_JDOS_M_spin[i] * norm_factor[q + 1] / sum_sum_JDOS_M_spin;

        int hits = 0;
        for (int i = 0; i < NE * line_size_Npos.at(q); i++)
            if (JDOS_M_spin[q][i] > 0)
                hits++;

        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        now = time(0);
        t = ctime(&now); t.pop_back();

        console_output = t + " | q: " + to_string(q) + "/" + to_string(q_max) + " | q_time: " + to_string(q_time) + "s | E: " + to_string(hits) + " | q_time/E: " + to_string(q_time / hits) + "s";
        string data_line = to_string(q) + " " + to_string(q_max) + " " + to_string(q_time) + " " + to_string(hits) + " " + to_string(q_time / hits);
        
        console_log.push_back(console_output);
        data.push_back(data_line);

        cout << console_output << endl;
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
    fs::create_directories(SAVE_DIR(lattice, L, (int) log10(REP)));

    std::ofstream file1((string) SAVE_DIR(lattice, L, (int) log10(REP)) + save_file + ".txt");
    for (int i = 0; i < NE; i++) 
    {
        for (int j = 0; j < NM; j++) 
            file1 << JDOS[i * NM + j] << " ";
        file1 << "\n";
    }
    file1.close();

    std::ofstream file2((string) SAVE_DIR(lattice, L, (int) log10(REP)) + save_file + "_data.txt");
    file2 << "q q_max q_time hits q_time/hits \n"; 
    for (int i = 0; i < data.size(); i++)
        file2 << data[i] << "\n";
    file2 << runtime << "\n";
    file2.close();

    std::ofstream file3((string) SAVE_DIR(lattice, L, (int) log10(REP)) + save_file + "_console_logs.txt");
    for (int i = 0; i < console_log.size(); i++)
        file3 << console_log.at(i) << "\n";
    file3.close();
    
    // Deallocate arrays

    for (int i = 0; i < line_size_Npos.size(); i++)
        delete[] Npos[i];

    for (int i = 0; i < NM; i++)
        delete[] JDOS_M_spin[i];

    delete[] JDOS_M_spin;
    delete[] JDOS, delete[] hist, delete[] hist_E_selected;
    delete[] new_spins_vector, delete[] spins_vector;
    delete[] NN_table, delete[] norm_factor, delete[] Npos;

    return 0;
}

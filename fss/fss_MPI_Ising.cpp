//
// MPI implementation of the Flat Scan Sampling for the Ising 1/2 Model 
// João Inácio, Apr. 17th, 2021
//
// This version is parallelized with MPI
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
#include <mpi.h>

#include "Fss_Functions.h"
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
#define SAVE_DIR(lattice, L, log_REP, workers)    "./data/mc/n_cores" + to_string(workers) + "/" + lattice + "/L" + to_string(L) + "/" + to_string(log_REP) + "/"


int main(int argc, char **argv)
{
    // Root, rank and size

    int root = 0;
    int rank;
    int size;

    // Initialize MPI 

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int workers = size - 1;

    // Set the seed for xoshiro256++

    uint64_t seed = SEED;
    if (seed == 0)
    {
        srand((unsigned) time(NULL));
        seed = rand();
    }

    splitmix64_seed(seed * rank);
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

    int q_max = (NM + 1) / 2 - 2;
    if (NM % 2 == 0)
        q_max = NM / 2 - 3;
    
    int run;
    int skip;
    ll REP;
    ll shuffle_REP;

    switch (argc)
    {
        case 2:
            REP = pow(10, atoi(argv[1]));
            run = 0;
            skip = N_atm;
            shuffle_REP = REP;
            break;
        case 3:
            REP = pow(10, atoi(argv[1]));
            run = atoi(argv[2]);
            skip = N_atm;
            shuffle_REP = REP;
            break;
        case 4:
            REP = pow(10, atoi(argv[1]));
            shuffle_REP = pow(10, atoi(argv[2]));
            run = atoi(argv[3]);
            skip = N_atm;
            break;
        case 5:
            REP = pow(10, atoi(argv[1]));
            shuffle_REP = pow(10, atoi(argv[2]));
            skip = atoi(argv[3]);
            run = atoi(argv[4]);
            break;
        default:
            if (rank == root)
                cout << "No parameters selected, reseting do default." << endl;
            run = 0;
            skip = N_atm;
            REP = pow(10, 4);
            shuffle_REP = REP;
            break;
    }

    ll REP_worker = REP / workers;

    string NN_table_file_name = "./neighbour_tables/neighbour_table_" + to_string(dim) + "D_" + 
    lattice + "_" + to_string(NN) + "NN_L" + to_string(L) + ".txt";
    string norm_factor_file_name = "./coefficients/coefficients_" + to_string(N_atm) + "d2.txt";
    string save_file = to_string(run) + "_JDOS_FSS_Ising_" + to_string(dim) + "D_" + lattice + "_L" + to_string(L) + 
    "_REP_1E" + to_string((int) log10(REP)) + "_skip_" + to_string(skip) + "_shuffle_" + to_string((int) log10(shuffle_REP));

    // Initialize vectors and read files

    int *spins_vector = new int[N_atm];
    int *new_spins_vector = new int[N_atm];
    int *NN_table = new int[N_atm * NN];
    ld *norm_factor = new ld[NM];
    
    read_NN_talbe(NN_table_file_name, NN_table);
    read_norm_factor(norm_factor_file_name, norm_factor);

    ld *JDOS = new ld[NE * NM];
    ld *JDOS_root = new ld[NE * NM];
    ll *hist = new ll[NE];
    ll *hist_E_selected = new ll[NE];

    for (int i = 0; i < NE * NM; i++)
        JDOS[i] = 0;
    JDOS[0] = 1;

    int q;
    int hits;
    int hits_root;

    // Start measuring time

    vector<string> console_log;
    vector<string> data;
    time_t now;
    string t;
    
    if (rank == root)
    {
        now = time(0);
        t = ctime(&now); t.pop_back();

        string console_output = "run: " + to_string(run) + " | L: " + to_string(L) + " | REP: " + to_string(REP) + " | skip: " + 
        to_string(skip) + " | dim: " + to_string(dim) + "D | lattie: " + lattice;
        console_log.push_back(console_output);
        string console_output_2 = "shuffle: " + to_string(shuffle_REP) + " | workers: " + to_string(workers) + " | REP/worker: " + 
        to_string(REP_worker);
        console_log.push_back(console_output_2);

        cout << endl;
        cout << console_output << endl;
        cout << console_output_2 << endl;
        cout << "Starting time: " << t << endl << endl;

        auto method_start = std::chrono::steady_clock::now();

        // Flat Scan Sampling
        // Scan and Compute JDOS at q = 1

        for (int flip_idx = 0; flip_idx < N_atm; flip_idx++)
        {
            int E_tmp1 = 0;
            for (int a = 0; a < NN; a++)
                E_tmp1 += - 1;
            
            int E_tmp2 = 0;
            for (int a = 0; a < NN; a++)
                E_tmp2 += 1;
            
            int E_tmp3 = - max_E - E_tmp1 + E_tmp2;
            int idx_E_tmp3 = energies[E_tmp3];

            JDOS[idx_E_tmp3 * NM + 1] += JDOS[0];
        }

        int sum_JDOS = 0;
        for (int i = 0; i < NE; i++)
            if (JDOS[i * NM + 1] > 0)
                sum_JDOS += JDOS[i * NM + 1];

        for (int i = 0; i < NE; i++)
            JDOS[i * NM + 1] = JDOS[i * NM + 1] * norm_factor[1] / sum_JDOS;

        console_output = t + " | q: " + to_string(0) + "/" + to_string(q_max);
        console_log.push_back(console_output);

        cout << console_output << endl;

        // Main Loop manager

        q = 1;
        MPI_Bcast(&q, 1, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(JDOS, NE * NM, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);

        for (q = 1; q <= q_max; q++)
        {
            auto q_start = std::chrono::steady_clock::now();

            for (int i = 0; i < NE * NM; i++)
                JDOS_root[i] = 0;
            hits_root = 0;
            hits = 0;
            
            MPI_Reduce(JDOS, JDOS_root, NE * NM, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            MPI_Reduce(&hits, &hits_root, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

            // Normalize JDOS and output to console

            hits = hits_root / workers;

            ld sum_JDOS = 0;
            for (int i = 0; i < NE; i++)
                if (JDOS_root[i * NM + q + 1] > 0)
                    sum_JDOS += JDOS_root[i * NM + q + 1];

            for (int i = 0; i < NE; i++)
                JDOS[i * NM + q + 1] = JDOS_root[i * NM + q + 1] * norm_factor[q + 1] / sum_JDOS;

            q++;
            if (q != q_max + 1)
            {
                MPI_Bcast(&q, 1, MPI_INT, root, MPI_COMM_WORLD);
                MPI_Bcast(JDOS, NE * NM, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);
            }
            q--;
            
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
        fs::create_directories(SAVE_DIR(lattice, L, (int) log10(REP), workers));

        std::ofstream file1((string) SAVE_DIR(lattice, L, (int) log10(REP), workers) + save_file + ".txt");
        for (int i = 0; i < NE; i++) 
        {
            for (int j = 0; j < NM; j++) 
                file1 << JDOS[i * NM + j] << " ";
            file1 << "\n";
        }
        file1.close();

        std::ofstream file2((string) SAVE_DIR(lattice, L, (int) log10(REP), workers) + save_file + "_data.txt");
        file2 << "q q_max q_time hits q_time/hits \n"; 
        for (int i = 0; i < data.size(); i++)
            file2 << data[i] << "\n";
        file2 << runtime << "\n";
        file2.close();

        std::ofstream file3((string) SAVE_DIR(lattice, L, (int) log10(REP), workers) + save_file + "_console_logs.txt");
        for (int i = 0; i < console_log.size(); i++)
            file3 << console_log.at(i) << "\n";
        file3.close();

    }
    else
    {
        // Flat Scan Sampling
        // Main Loop worker

        do
        {
            MPI_Bcast(&q, 1, MPI_INT, root, MPI_COMM_WORLD);
            MPI_Bcast(JDOS, NE * NM, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);

            for (int i = 0; i < NE; i++)
            {
                hist[i] = 0; 
                hist_E_selected[i] = 0;
            }

            // Random config at q

            for (int i = 0; i < N_atm; i++)
                spins_vector[i] = 1;
            int E_config =  - max_E;
            int idx_E_config = energies[E_config];

            array<vector<int>, 2> flip_list;
            for (int i = 0; i < N_atm; i++)
                flip_list[0].push_back(i);

            for (int idx = 1; idx <= q; idx++)
            {
                int idx_tmp = rand_xoshiro256pp() % flip_list[0].size();
                int flipped_idx = flip_list[0].at(idx_tmp);
                spins_vector[flipped_idx] = - 1;

                flip_list[1].push_back(flipped_idx);
                flip_list[0].erase(flip_list[0].begin() + idx_tmp);
                
                int delta_E = 0;
                for (int a = 0; a < NN; a++)
                    delta_E += - spins_vector[flipped_idx] * spins_vector[NN_table[flipped_idx * NN + a]];
                
                E_config += 2 * delta_E;
            }
            idx_E_config = energies[E_config];
            
            // Shuffle_configuration

            shuffle(JDOS, shuffle_REP, flip_list, spins_vector, q, N_atm, NN, NM, NN_table, E_config, idx_E_config, energies);

            // Update Histograms

            hist[idx_E_config]++;
            hist_E_selected[idx_E_config]++;

            // Scan the first config

            for (int flip_idx = 0; flip_idx < flip_list[0].size(); flip_idx++)
            {
                int delta_E = 0;
                for (int a = 0; a < NN; a++)
                    delta_E += spins_vector[flip_list[0].at(flip_idx)] * spins_vector[NN_table[flip_list[0].at(flip_idx) * NN + a]];

                int E_tmp = E_config + 2 * delta_E;
                int idx_E_tmp = energies[E_tmp];

                JDOS[idx_E_tmp * NM + q + 1] += JDOS[idx_E_config * NM + q] / REP_worker;
            }

            ll k = 1;
            bool accepted = false;

            // Where the magic happens

            while (min_hist(hist_E_selected, NE) < REP_worker)
            {
                // Get a new random condig at magnetization q

                if (!accepted)
                    for (int i = 0; i < N_atm; i++)
                        new_spins_vector[i] =  spins_vector[i];
                int new_E_config = 0;
                int new_idx_E_config = 0;

                // Flip a positive spin to a negative

                int idx_tmp1 = rand_xoshiro256pp() % flip_list[0].size();
                int flipped_idx1 = flip_list[0].at(idx_tmp1);
                new_spins_vector[flipped_idx1] = - 1;

                flip_list[1].push_back(flipped_idx1);
                flip_list[0].erase(flip_list[0].begin() + idx_tmp1);

                int delta_E = 0;
                for (int a = 0; a < NN; a++)
                    delta_E += - new_spins_vector[flipped_idx1] * new_spins_vector[NN_table[flipped_idx1 * NN + a]];
                new_E_config = E_config + 2 * delta_E;

                // Flip a negative spin to a positive

                int idx_tmp2 = rand_xoshiro256pp() % flip_list[1].size();
                int flipped_idx2 = flip_list[1].at(idx_tmp2);
                new_spins_vector[flipped_idx2] = 1;

                flip_list[0].push_back(flipped_idx2);
                flip_list[1].erase(flip_list[1].begin() + idx_tmp2);

                delta_E = 0;
                for (int a = 0; a < NN; a++)
                    delta_E += - new_spins_vector[flipped_idx2] * new_spins_vector[NN_table[flipped_idx2 * NN + a]];
                new_E_config = new_E_config + 2 * delta_E;          

                // Wang Landau criteria

                new_idx_E_config = energies[new_E_config];
                ld ratio = JDOS[idx_E_config * NM + q] / JDOS[new_idx_E_config * NM + q];

                if (ratio >= 1 || ((ld) rand_xoshiro256pp() / (ld) UINT64_MAX) < ratio || hist_E_selected[new_idx_E_config] == 0)
                {
                    for (int i = 0; i < N_atm; i++)
                        spins_vector[i] = new_spins_vector[i];

                    E_config = new_E_config;
                    idx_E_config = new_idx_E_config;
                    accepted = true;
                }
                else
                {
                    if (flipped_idx1 != flipped_idx2)
                    {
                        flip_list[0].pop_back();
                        flip_list[0].push_back(flipped_idx1);

                        flip_list[1].pop_back();
                        flip_list[1].push_back(flipped_idx2);
                    }
                    accepted = false;
                }

                hist[idx_E_config]++;

                // Scan configuration

                if (hist_E_selected[idx_E_config] < REP_worker && k % skip == 0 || hist_E_selected[new_idx_E_config] == 0)
                {
                    for (int flip_idx = 0; flip_idx < flip_list[0].size(); flip_idx++)
                    {
                        int delta_E = 0;
                        for (int a = 0; a < NN; a++)
                            delta_E += spins_vector[flip_list[0].at(flip_idx)] * spins_vector[NN_table[flip_list[0].at(flip_idx) * NN + a]];

                        int E_tmp = E_config + 2 * delta_E;
                        int idx_E_tmp = energies[E_tmp];

                        JDOS[idx_E_tmp * NM + q + 1] += JDOS[idx_E_config * NM + q] / REP_worker;
                    }

                    hist_E_selected[idx_E_config]++;
                }

                k++;
            }

            // Compute hits and send results

            hits = 0;
            for (int i = 0; i < NE; i++)
                if (JDOS[i * NM + q] > 0)
                    hits++;

            MPI_Reduce(JDOS, JDOS_root, NE * NM, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            MPI_Reduce(&hits, &hits_root, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
        } while (q != q_max);
    }
    
    delete[] JDOS, hist, hist_E_selected, JDOS_root;
    delete[] new_spins_vector, spins_vector;
    delete[] NN_table, norm_factor;

    MPI_Finalize();

    return 0;
}


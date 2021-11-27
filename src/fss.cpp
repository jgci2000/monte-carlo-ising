
#include "fss.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <stdint.h>

FSS::FSS(RNG &rng, System &ising_lattice) {
    this->set_rng(rng);
    this->set_lattice(ising_lattice);
}

FSS::FSS(RNG &rng) {
    this->set_rng(rng);
}

FSS::FSS(long long REP, int skip, RNG &rng) {
    this->set_rng(rng);
    this->set_params(REP, skip);
}

FSS::FSS(long long REP, int skip, RNG &rng, System &ising_lattice) {
    this->set_rng(rng);
    this->set_params(REP, skip);
    this->set_lattice(ising_lattice);
}

FSS::~FSS() {
    delete[] this->hist;
    delete[] this->time_iter;
    delete[] this->steps_iter;
}

void FSS::set_params(long long REP, int skip) {
    this->REP = REP;
    this->skip = skip;

    this->run_time = 0;
    
    this->added_params = true;
}

void FSS::set_lattice(System &ising_lattice) {
    this->ising_lattice = &(ising_lattice);
    this->ising_lattice->init_spins_random(*(this->rng));

    this->hist = new long long[this->ising_lattice->NE];

    for (int i = 0; i < this->ising_lattice->NE; ++i) {
        this->hist[i] = 0;
    }

    this->idx_M_max = (this->ising_lattice->NM + 1) / 2 - 2;
    if (this->ising_lattice->NM % 2 == 0)
        this->idx_M_max = this->ising_lattice->NM / 2 - 3;

    this->time_iter = new double[this->idx_M_max];
    this->steps_iter = new long long[this->idx_M_max];

    this->added_lattice = true;
}

void FSS::set_rng(RNG &rng) {
    this->rng = &(rng);

    this->added_rng = true;
}

void FSS::simulate(int run, bool verbose) {
    if (!(this->added_lattice && this->added_params && this->added_rng)) {
        std::printf(" -- Error: forgot to add the simulation parameters, rng or lattice -- ");
    }

    std::printf("Initiating Flat Scan Sampling Simulation; run: %d \n", run);
    time_t now = time(0);
    std::string t = ctime(&now); t.pop_back();
    std::printf("    Time: %s \n", t.c_str());
    if (verbose) {    
        std::printf("    System:  L: %d | Sz: %d | N_atm: %d | lattice: %s | NN: %d \n",
            this->ising_lattice->L, 
            this->ising_lattice->Sz,
            this->ising_lattice->N_atm,
            this->ising_lattice->lattice.c_str(), 
            this->ising_lattice->NN);
        std::printf("    Simulation Parameters: REP: %f | skip: %f \n",
            this->REP, 
            this->skip);   
        std::printf("\n");
    }

    // init variables
    int *new_spins_vector = new int[this->ising_lattice->N_atm];

    auto runtime_start = std::chrono::steady_clock::now();

    this->ising_lattice->JDOS[0] = 1;

    for (int flip_idx = 0; flip_idx < this->ising_lattice->N_atm; flip_idx++)
    {
        int E_tmp1 = 0;
        for (int a = 0; a < this->ising_lattice->NN; a++)
            E_tmp1 += - 1;
        
        int E_tmp2 = 0;
        for (int a = 0; a < this->ising_lattice->NN; a++)
            E_tmp2 += 1;
        
        int E_tmp3 = - this->ising_lattice->max_E - E_tmp1 + E_tmp2;
        int idx_E_tmp3 = this->ising_lattice->energies[E_tmp3];

        this->ising_lattice->JDOS[idx_E_tmp3 * this->ising_lattice->NM + 1] += this->ising_lattice->JDOS[0];
    }

    int sum_JDOS = 0;
    for (int i = 0; i < this->ising_lattice->NE; i++)
        if (this->ising_lattice->JDOS[i * this->ising_lattice->NM + 1] > 0)
            sum_JDOS += this->ising_lattice->JDOS[i * this->ising_lattice->NM + 1];

    for (int i = 0; i < this->ising_lattice->NE; i++)
        this->ising_lattice->JDOS[i * this->ising_lattice->NM + 1] = this->ising_lattice->JDOS[i * this->ising_lattice->NM + 1] * exp(this->ising_lattice->norm_factor[1]) / sum_JDOS;

    // LOOP

    for (int q = 1; q <= this->idx_M_max; ++q) {
        auto q_start = std::chrono::steady_clock::now();

        for (int i = 0; i < this->ising_lattice->NE; ++i) {
            this->hist[i] = 0; 
        }

        // Random config at q

        for (int i = 0; i < this->ising_lattice->N_atm; i++)
            this->ising_lattice->spins_vector[i] = 1;
        int E_config =  - this->ising_lattice->max_E;
        int idx_E_config = this->ising_lattice->energies[E_config];

        std::array<std::vector<int>, 2> flip_list;
        for (int i = 0; i < this->ising_lattice->N_atm; i++)
            flip_list[0].push_back(i);

        for (int idx = 1; idx <= q; idx++)
        {
            int idx_tmp = this->rng->rand() % flip_list[0].size();
            int flipped_idx = flip_list[0].at(idx_tmp);
            this->ising_lattice->spins_vector[flipped_idx] = - 1;

            flip_list[1].push_back(flipped_idx);
            flip_list[0].erase(flip_list[0].begin() + idx_tmp);
            
            int delta_E = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                delta_E += - this->ising_lattice->spins_vector[flipped_idx] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flipped_idx * this->ising_lattice->NN + a]];
            
            E_config += 2 * delta_E;
        }
        idx_E_config = this->ising_lattice->energies[E_config];

        // Update Histograms

        this->hist[idx_E_config]++;

        // Scan the first config

        for (int flip_idx = 0; flip_idx < flip_list[0].size(); flip_idx++)
        {
            int delta_E = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                delta_E += this->ising_lattice->spins_vector[flip_list[0].at(flip_idx)] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flip_list[0].at(flip_idx) * this->ising_lattice->NN + a]];

            int E_tmp = E_config + 2 * delta_E;
            int idx_E_tmp = this->ising_lattice->energies[E_tmp];

            this->ising_lattice->JDOS[idx_E_tmp * this->ising_lattice->NM + q + 1] += this->ising_lattice->JDOS[idx_E_config * this->ising_lattice->NM + q] / this->REP;
        }

        long long k = 1;
        bool accepted = false;

        // Where the magic happens

        while (min_hist() < this->REP)
        {
            // Get a new random condig at magnetization q

            if (!accepted)
                for (int i = 0; i < this->ising_lattice->N_atm; i++)
                    new_spins_vector[i] =  this->ising_lattice->spins_vector[i];
            int new_E_config = 0;
            int new_idx_E_config = 0;

            // Flip a positive spin to a negative

            int idx_tmp1 = this->rng->rand() % flip_list[0].size();
            int flipped_idx1 = flip_list[0].at(idx_tmp1);
            new_spins_vector[flipped_idx1] = - 1;

            flip_list[1].push_back(flipped_idx1);
            flip_list[0].erase(flip_list[0].begin() + idx_tmp1);

            int delta_E = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                delta_E += - new_spins_vector[flipped_idx1] * new_spins_vector[this->ising_lattice->NN_table[flipped_idx1 * this->ising_lattice->NN + a]];
            new_E_config = E_config + 2 * delta_E;

            // Flip a negative spin to a positive

            int idx_tmp2 = this->rng->rand() % flip_list[1].size();
            int flipped_idx2 = flip_list[1].at(idx_tmp2);
            new_spins_vector[flipped_idx2] = 1;

            flip_list[0].push_back(flipped_idx2);
            flip_list[1].erase(flip_list[1].begin() + idx_tmp2);

            delta_E = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                delta_E += - new_spins_vector[flipped_idx2] * new_spins_vector[this->ising_lattice->NN_table[flipped_idx2 * this->ising_lattice->NN + a]];
            new_E_config = new_E_config + 2 * delta_E;          

            // Wang Landau criteria

            new_idx_E_config = this->ising_lattice->energies[new_E_config];
            long double ratio = this->ising_lattice->JDOS[idx_E_config * this->ising_lattice->NM + q] / this->ising_lattice->JDOS[new_idx_E_config * this->ising_lattice->NM + q];

            if (ratio >= 1 || this->rng->rand_uniform() < ratio || this->hist[new_idx_E_config] == 0)
            {
                for (int i = 0; i < this->ising_lattice->N_atm; i++)
                    this->ising_lattice->spins_vector[i] = new_spins_vector[i];

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

            this->hist[idx_E_config]++;

            // Scan configuration

            if (this->hist[idx_E_config] < REP && k % skip == 0 || this->hist[new_idx_E_config] == 0)
            {
                for (int flip_idx = 0; flip_idx < flip_list[0].size(); flip_idx++)
                {
                    int delta_E = 0;
                    for (int a = 0; a < this->ising_lattice->NN; a++)
                        delta_E += this->ising_lattice->spins_vector[flip_list[0].at(flip_idx)] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flip_list[0].at(flip_idx) * this->ising_lattice->NN + a]];

                    int E_tmp = E_config + 2 * delta_E;
                    int idx_E_tmp = this->ising_lattice->energies[E_tmp];

                    this->ising_lattice->JDOS[idx_E_tmp * this->ising_lattice->NM + q + 1] += this->ising_lattice->JDOS[idx_E_config * this->ising_lattice->NM + q] / this->REP;
                }

                this->hist[idx_E_config]++;
            }

            k++;
        }

        // Normalize JDOS and output to console

        long double sum_JDOS = 0;
        for (int i = 0; i < this->ising_lattice->NE; i++)
            if (this->ising_lattice->JDOS[i * this->ising_lattice->NM + q + 1] > 0)
                sum_JDOS += this->ising_lattice->JDOS[i * this->ising_lattice->NM + q + 1];

        for (int i = 0; i < this->ising_lattice->NE; i++)
            this->ising_lattice->JDOS[i * this->ising_lattice->NM + q + 1] = this->ising_lattice->JDOS[i * this->ising_lattice->NM + q + 1] * exp(this->ising_lattice->norm_factor[q + 1]) / sum_JDOS;

        int hits = 0;
        for (int i = 0; i < this->ising_lattice->NE; i++)
            if (this->ising_lattice->JDOS[i * this->ising_lattice->NM + q] > 0)
                hits++;

        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        now = time(0);
        t = ctime(&now); t.pop_back();

        std::printf("done q=%d \n", q);
    }

    // LOOP

    auto runtime_end = std::chrono::steady_clock::now();
    this->run_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (runtime_end - runtime_start).count()) * pow(10, -6);

    if (verbose) {
        std::printf("\n");
    }
    std::printf("    Run time: %fs \n", this->run_time);
    now = time(0);
    t = ctime(&now); t.pop_back();
    std::printf("    Time: %s \n", t.c_str());
    std::printf("Finished Flat Scan Sampling Simulation; run : %d \n", run);

    delete[] new_spins_vector;
}

void FSS::write_to_file(std::string name, bool debug) {
    int status = system("mkdir results");

    std::ofstream file1("results/" + name);
    if (file1.is_open()) {
        for (int i = 0; i < this->ising_lattice->NE; ++i) 
        {
            for (int j = 0; j < this->ising_lattice->NM; ++j) {
                file1 << this->ising_lattice->JDOS[i * this->ising_lattice->NM + j] << " ";
            }
            file1 << "\n";
        }
        file1.close();
    } else {
        std:printf(" -- Error: can not open save file, please check you directory -- \n");
    }

    if (debug) {
        int status = system("mkdir results/debug");

        std::ofstream file2("results/debug/" + name);
        if (file2.is_open()) {
            file2 << this->run_time << "\n";
            for (int i = 0; i < this->n_f_vals; i++) {
                file2 << this->steps_iter[i] << " " << this->time_iter[i] << "\n";
            }
            file2.close();
        } else {
            std::printf(" -- Error: can not open debug file, please check you directory -- \n");
        }
    }
}

void FSS::print_JDOS() {
    std::printf("\nJDOS: \n");
    for (int i = 0; i < this->ising_lattice->NE; ++i) 
    {
        for (int j = 0; j < this->ising_lattice->NM; ++j) {
            std::cout << this->ising_lattice->JDOS[i * this->ising_lattice->NM + j] << " ";
        }
        std::cout << std::endl;
    } 
}

void FSS::normalize_JDOS() {
    
}

long long FSS::min_hist() {
    long long min = __LONG_LONG_MAX__;
    for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; i++)
        if (hist[i] != 0 && this->hist[i] < min)
            min = this->hist[i];
    return min;
}


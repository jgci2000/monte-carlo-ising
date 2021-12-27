/**
 * Implementation by João Inácio (j.g.c.inacio@fys.uio.no).
 * Dec. 2021
 */

#include "wl.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdlib.h>

/**
 * Wang Landau class. 
 * @param rng (RNG): Reference to an RNG object.
 * @param ising_lattice (System): reference to a System object. 
 */
WL::WL(RNG &rng, System &ising_lattice) {
    this->set_rng(rng);
    this->set_lattice(ising_lattice);
}

/**
 * Wang Landau class. 
 * @param rng (RNG): Reference to an RNG object.
 */
WL::WL(RNG &rng) {
    this->set_rng(rng);
}

/**
 * Wang Landau class. 
 * @param f_init (double): Initial value for the inital modification factor. Usually set to exp(1).
 * @param f_final (double): Final value for the modification factor. Should be very close to 1, e.g. 1 + 10^(-8).
 * @param flatness (double): Flatness criteria. Value between 0 and 1. Usually 0.9.
 * @param rng (RNG): Reference to an RNG object.
 */
WL::WL(double f_init, double f_final, double flatness, RNG &rng) {
    this->set_rng(rng);
    this->set_params(f_init, f_final, flatness);
}

/**
 * Wang Landau class. 
 * @param f_init (double): Initial value for the inital modification factor. Usually set to exp(1).
 * @param f_final (double): Final value for the modification factor. Should be very close to 1, e.g. 1 + 10^(-8).
 * @param flatness (double): Flatness criteria. Value between 0 and 1. Usually 0.9.
 * @param rng (RNG): Reference to an RNG object.
 * @param ising_lattice (System): reference to a System object. 
 */
WL::WL(double f_init, double f_final, double flatness, RNG &rng, System &ising_lattice) {
    this->set_rng(rng);
    this->set_params(f_init, f_final, flatness);
    this->set_lattice(ising_lattice);
}

WL::~WL() {
    delete[] this->ln_JDOS;
    delete[] this->hist;
    delete[] this->time_iter;
    delete[] this->steps_iter;
}

/**
 * Sets parameters for Wang Landau class. 
 * @param f_init (double): Initial value for the inital modification factor. Usually set to exp(1).
 * @param f_final (double): Final value for the modification factor. Should be very close to 1, e.g. 1 + 10^(-8).
 * @param flatness (double): Flatness criteria. Value between 0 and 1. Usually 0.9.
 */
void WL::set_params(double f_init, double f_final, double flatness) {
    this->f_init = f_init;
    this->f_final = f_final;
    this->flatness = flatness;

    this->n_f_vals = 0;
    while (f_init > this->f_final) {
        n_f_vals++;
        f_init = sqrt(f_init);
    }

    this->time_iter = new double[this->n_f_vals];
    this->steps_iter = new unsigned long long[this->n_f_vals];

    this->run_time = 0;
    
    this->added_params = true;
}

/**
 * Sets Ising lattice for Wang Landau class.
 * @param ising_lattice (System): Reference to System object.
 */
void WL::set_lattice(System &ising_lattice) {
    this->ising_lattice = &(ising_lattice);
    this->ising_lattice->init_spins_random(*(this->rng));

    this->ln_JDOS = new long double[this->ising_lattice->NE * this->ising_lattice->NM];
    this->hist = new unsigned long long[this->ising_lattice->NE * this->ising_lattice->NM];

    for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; ++i) {
        this->ln_JDOS[i] = 0;
        this->hist[i] = 0;
    }

    this->added_lattice = true;
}

/**
 * Sets RNG for the Wang Landau class.
 * @param rng (RNG): Reference to RNG object.
 */
void WL::set_rng(RNG &rng) {
    this->rng = &(rng);

    this->added_rng = true;
}

/**
 * Estimates the Joint Density of States for with the Wang Landau
 * algorithm, with respect to the parameters and lattice passes beforehand.
 * @param steps (unsigned long): Number of steps taken before checking for histogram flatness.
 * @param run (int): Should be passed 0, if only doing one run. Used to keep track of how many times 
 * the algorithm has run in a given script.
 * @param verbose (bool): If true, output information about the simulation each iteration.
 */
void WL::simulate(unsigned long long steps, int run, bool verbose) {
    if (!(this->added_lattice && this->added_params && this->added_rng)) {
        printf(" -- Error: forgot to add the simulation parameters, rng or lattice -- ");
    }

    printf("Initiating Wang-Landau Simulation; run: %d \n", run);
    time_t now = time(0);
    std::string t = ctime(&now); t.pop_back();
    printf("    Time: %s \n", t.c_str());
    if (verbose) {    
        printf("    System:  L: %d | Sz: %d | N_atm: %d | lattice: %s | NN: %d \n",
            this->ising_lattice->L, 
            this->ising_lattice->Sz,
            this->ising_lattice->N_atm,
            this->ising_lattice->lattice.c_str(), 
            this->ising_lattice->NN);
        printf("    Simulation Parameters: f_init: %f | f_final: %.11f | flatness: %.2f | steps %ld \n",
            this->f_init, 
            this->f_final, 
            this->flatness, 
            steps);   
        printf("\n");
    }

    bool take_time = true;

    double f = this->f_init;

    int idx_E_config = this->ising_lattice->energies[this->ising_lattice->E_config];
    int idx_M_config = this->ising_lattice->magnetizations[this->ising_lattice->M_config];

    int **flip_to = new int*[this->ising_lattice->N_atm];
    int *is_at = new int[this->ising_lattice->N_atm];
    for (int i = 0; i < this->ising_lattice->N_atm; ++i) {
        int c = 0;
        flip_to[i] = new int[this->ising_lattice->Sz - 1];
        for (int j = 0; j < this->ising_lattice->Sz; ++j) {
            if (this->ising_lattice->spins_vector[i] != this->ising_lattice->spins_values[j]) {
                flip_to[i][c] = j;
                c++;
            }
            else {
                is_at[i] = j;
            }
        }
    }

    int flip_idx, flip_spin_to_idx, tmp_idx;
    int delta_E, delta_M;
    int new_E_config, new_M_config;
    int new_idx_E_config, new_idx_M_config;

    long double ratio;
    
    long long mc_sweep = 0;
    int counter = 0;

    auto loop_start = std::chrono::steady_clock::now();
    auto runtime_start = std::chrono::steady_clock::now();

    while (f > this->f_final) {
        if (take_time) {
            loop_start = std::chrono::steady_clock::now();
            take_time = false;
        }

        for (int idx = 0; idx < this->ising_lattice->N_atm * steps; ++idx) {
            flip_idx = this->rng->rand() % this->ising_lattice->N_atm;
            tmp_idx = this->rng->rand() % (this->ising_lattice->Sz - 1);
            flip_spin_to_idx = flip_to[flip_idx][tmp_idx];

            delta_E = this->ising_lattice->energy_flip(flip_idx, flip_spin_to_idx);
            delta_M = this->ising_lattice->magnetization_flip(flip_idx, flip_spin_to_idx);

            new_E_config = this->ising_lattice->E_config + delta_E;
            new_M_config  = this->ising_lattice->M_config + delta_M;
            new_idx_E_config = this->ising_lattice->energies[new_E_config];
            new_idx_M_config = this->ising_lattice->magnetizations[new_M_config];

            ratio = exp(ln_JDOS[idx_E_config * this->ising_lattice->NM + idx_M_config] - ln_JDOS[new_idx_E_config * this->ising_lattice->NM + new_idx_M_config]);

            if (ratio >= 1 || this->rng->rand_uniform() < ratio) {
                this->ising_lattice->spins_vector[flip_idx] = this->ising_lattice->spins_values[flip_spin_to_idx];

                this->ising_lattice->E_config = new_E_config;
                idx_E_config = new_idx_E_config;
                this->ising_lattice->M_config = new_M_config;
                idx_M_config = new_idx_M_config;

                flip_to[flip_idx][tmp_idx] = is_at[flip_idx];
                is_at[flip_idx] = flip_spin_to_idx;
            }

            this->hist[idx_E_config * this->ising_lattice->NM + idx_M_config]++;
            this->ln_JDOS[idx_E_config * this->ising_lattice->NM + idx_M_config] += log(f);
        }
        mc_sweep += steps;

        if (this->flat_hist()) {
            auto loop_end = std::chrono::steady_clock::now();
            double loop_dur = (double) (std::chrono::duration_cast<std::chrono::microseconds> (loop_end - loop_start).count()) * pow(10, -6);

            counter++;
            take_time = true;

            if (verbose) {
                now = time(0);
                t = ctime(&now); t.pop_back();
                printf("%s | run: %d | f: %d/%d | time: %fs | steps: %ld \n", 
                    t.c_str(),
                    run, 
                    counter, 
                    this->n_f_vals, 
                    loop_dur, 
                    mc_sweep);
            }
            
            this->time_iter[counter - 1] = loop_dur;
            this->steps_iter[counter - 1] = mc_sweep;

            f = sqrt(f);
            mc_sweep = 0;
            for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; ++i) {
                this->hist[i] = 0;
            }
        }
    }
    this->normalize_JDOS();

    auto runtime_end = std::chrono::steady_clock::now();
    this->run_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (runtime_end - runtime_start).count()) * pow(10, -6);

    if (verbose) {
        printf("\n");
    }
    printf("    Run time: %fs \n", this->run_time);
    now = time(0);
    t = ctime(&now); t.pop_back();
    printf("    Time: %s \n", t.c_str());
    printf("Finished Wang-Landau Simulation; run : %d \n", run);

    for (int i = 0; i < this->ising_lattice->N_atm; ++i) {
        delete[] flip_to[i];
    }
    delete[] flip_to;
    delete[] is_at;
}


/**
 * Writes estimated Joint Density of States to file.
 * @param name (string): Name of the file with extension.
 * @param path (string): Path of the directory where the file should be written to. 
 * Relative path from the execution directory, of the form "folder1/folder2/.../".
 * @param debug (bool): If true outputs additional information about the simulation.
 * First column is the number of steps and second is the time taken.
 */
void WL::write_to_file(std::string name, std::string path, bool debug) {
    std::filesystem::create_directories(path); 

    std::ofstream file1(path + name);
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
        std::filesystem::create_directories(path + "debug/"); 

        std::ofstream file2(path + "debug/" + name);
        if (file2.is_open()) {
            file2 << this->run_time << "\n";
            for (int i = 0; i < this->n_f_vals; i++) {
                file2 << this->steps_iter[i] << " " << this->time_iter[i] << "\n";
            }
            file2.close();
        } else {
            printf(" -- Error: can not open debug file, please check you directory -- \n");
        }
    }
}

/**
 * Prints the estimated Joint Density of States to the terminal.
 */
void WL::print_JDOS() {
    printf("\nJDOS: \n");
    for (int i = 0; i < this->ising_lattice->NE; ++i) 
    {
        for (int j = 0; j < this->ising_lattice->NM; ++j) {
            std::cout << this->ising_lattice->JDOS[i * this->ising_lattice->NM + j] << " ";
        }
        std::cout << std::endl;
    } 
}

void WL::normalize_JDOS() {
    for (int q = 0; q < this->ising_lattice->NM; ++q) {
        int first_idx;
        for (int i = 0; i < this->ising_lattice->NE; ++i) {
            if (this->ln_JDOS[i * this->ising_lattice->NM + q] > 0) {
                first_idx = i;
                break;
            }
        }

        long double temp = 0;
        for (int i = 0; i < this->ising_lattice->NE; ++i) {
            if (ln_JDOS[i * this->ising_lattice->NM + q] > 0) {
                temp += exp(this->ln_JDOS[i * this->ising_lattice->NM + q] - this->ln_JDOS[first_idx * this->ising_lattice->NM + q]);
            }
        }

        long double sum_ln_JDOS = ln_JDOS[first_idx * this->ising_lattice->NM + q] + log(temp);

        for (int i = 0; i < this->ising_lattice->NE; ++i) {
            if (this->ln_JDOS[i * this->ising_lattice->NM + q] > 0) {
                this->ising_lattice->JDOS[i * this->ising_lattice->NM + q] = exp(this->ln_JDOS[i * this->ising_lattice->NM + q] + this->ising_lattice->norm_factor[q] - sum_ln_JDOS);
            }
        }
    }
}

bool WL::flat_hist() {
    return this->min_hist() >= this->average_hist() * flatness ? true : false;
}

long double WL::average_hist() {
    long double sum = 0;
    int nnz = 0;
    for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; i++) {
        if (this->hist[i] != 0) {
            sum += this->hist[i];
            nnz++;
        }
    }
    return sum / nnz;
}

unsigned long long WL::min_hist() {
    unsigned long long min = __UINT64_MAX__;
    for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; i++)
        if (hist[i] != 0 && this->hist[i] < min)
            min = this->hist[i];
    return min;
}

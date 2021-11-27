
#include "wl.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <stdlib.h>

WL::WL(RNG &rng, System &ising_lattice) {
    this->set_rng(rng);
    this->set_lattice(ising_lattice);
}

WL::WL(RNG &rng) {
    this->set_rng(rng);
}

WL::WL(double f_init, double f_final, double flatness, RNG &rng) {
    this->set_rng(rng);
    this->set_params(f_init, f_final, flatness);
}

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
    this->steps_iter = new long long[this->n_f_vals];

    this->run_time = 0;
    
    this->added_params = true;
}

void WL::set_lattice(System &ising_lattice) {
    this->ising_lattice = &(ising_lattice);
    this->ising_lattice->init_spins_random(*(this->rng));

    this->ln_JDOS = new long double[this->ising_lattice->NE * this->ising_lattice->NM];
    this->hist = new long long[this->ising_lattice->NE * this->ising_lattice->NM];

    for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; ++i) {
        this->ln_JDOS[i] = 0;
        this->hist[i] = 0;
    }

    this->added_lattice = true;
}

void WL::set_rng(RNG &rng) {
    std::printf("here \n");
    this->rng = &(rng);
    std::printf("here \n");


    this->added_rng = true;
}

void WL::simulate(long long steps, int run, bool verbose) {
    if (!(this->added_lattice && this->added_params && this->added_rng)) {
        std::printf(" -- Error: forgot to add the simulation parameters, rng or lattice -- ");
    }

    std::printf("Initiating Wang-Landau Simulation; run: %d \n", run);
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
        std::printf("    Simulation Parameters: f_init: %f | f_final: %f | flatness: %f | steps %ld \n",
            this->f_init, 
            this->f_final, 
            this->flatness, 
            steps);   
        std::printf("\n");
    }

    bool take_time = true;

    double f = this->f_init;

    int idx_E_config = this->ising_lattice->energies[this->ising_lattice->E_config];
    int idx_M_config = this->ising_lattice->magnetizations[this->ising_lattice->M_config];

    int flip_idx;
    int delta_E;
    int new_E_config, new_M_config;
    int new_idx_E_config, new_idx_M_config;

    long double ratio;
    
    long long mc_sweep = 0;
    int counter = 1;

    auto loop_start = std::chrono::steady_clock::now();
    auto runtime_start = std::chrono::steady_clock::now();

    while (f > this->f_final) {
        if (take_time) {
            loop_start = std::chrono::steady_clock::now();
            take_time = false;
        }

        for (int idx = 0; idx < this->ising_lattice->N_atm * steps; ++idx) {
            flip_idx = this->rng->rand() % this->ising_lattice->N_atm;
            delta_E = this->ising_lattice->energy_flip(flip_idx);

            new_E_config = this->ising_lattice->E_config - 2 * delta_E;
            new_M_config  = this->ising_lattice->M_config - 2 * this->ising_lattice->spins_vector[flip_idx];
            new_idx_E_config = this->ising_lattice->energies[new_E_config];
            new_idx_M_config = this->ising_lattice->magnetizations[new_M_config];

            ratio = exp(ln_JDOS[idx_E_config * this->ising_lattice->NM + idx_M_config] - ln_JDOS[new_idx_E_config * this->ising_lattice->NM + new_idx_M_config]);

            if (ratio >= 1 || this->rng->rand_uniform() < ratio) {
                this->ising_lattice->spins_vector[flip_idx] =
                    - this->ising_lattice->spins_vector[flip_idx];

                this->ising_lattice->E_config = new_E_config;
                idx_E_config = new_idx_E_config;
                this->ising_lattice->M_config = new_M_config;
                idx_M_config = new_idx_M_config;
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
                std::printf("%s | run: %d | f: %d/%d | time: %fs | steps: %ld \n", 
                    t.c_str(),
                    run, 
                    counter, 
                    this->n_f_vals, 
                    loop_dur, 
                    mc_sweep);
            }
            
            this->time_iter[counter - 2] = loop_dur;
            this->steps_iter[counter - 2] = mc_sweep;

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
        std::printf("\n");
    }
    std::printf("    Run time: %fs \n", this->run_time);
    now = time(0);
    t = ctime(&now); t.pop_back();
    std::printf("    Time: %s \n", t.c_str());
    std::printf("Finished Wang-Landau Simulation; run : %d \n", run);
}

void WL::write_to_file(std::string name, bool debug) {
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

void WL::print_JDOS() {
    std::printf("\nJDOS: \n");
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

long long WL::min_hist() {
    long long min = __LONG_LONG_MAX__;
    for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; i++)
        if (hist[i] != 0 && this->hist[i] < min)
            min = this->hist[i];
    return min;
}

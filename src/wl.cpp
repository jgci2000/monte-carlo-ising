
#include "wl.h"

#include <chrono>

WL::WL(double f_init, double f_final, double flatness, RNG &rng) {
    this->init(f_init, f_final, flatness, rng);
}

WL::WL(double f_init, double f_final, double flatness, RNG &rng, System &ising_lattice) {
    this->init(f_init, f_final, flatness, rng);
    this->add_lattice(ising_lattice);
}

WL::~WL() {
    delete[] this->ln_JDOS;
    delete[] this->hist;
}

void WL::init(double f_init, double f_final, double flatness, RNG &rng) {
    this->f_init = f_init;
    this->f_final = f_final;
    this->flatness = flatness;

    this->rng = &(rng);

    this->n_f_vals = 0;
    while (f_init > this->f_final) {
        n_f_vals++;
        f_init = sqrt(f_init);
    }
}

void WL::add_lattice(System &ising_lattice) {
    this->ising_lattice = &(ising_lattice);
    this->ising_lattice->init_spins_random(*(this->rng));

    this->ln_JDOS = new long double[this->ising_lattice->NE * this->ising_lattice->NM];
    this->hist = new long long[this->ising_lattice->NE * this->ising_lattice->NM];

    for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; ++i) {
        this->ln_JDOS[i] = 0;
        this->hist[i] = 0;
    }
}

void WL::simulate(long long steps) {
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

    while (f > this->f_final) {
        if (mc_sweep == 0)
            loop_start = std::chrono::steady_clock::now();

        for (int idx = 0; idx < this->ising_lattice->N_atm; ++idx) {
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

        mc_sweep++;

        if (mc_sweep % steps == 0 && this->flat_hist()) {
            auto loop_end = std::chrono::steady_clock::now();
            double loop_dur = (double) (std::chrono::duration_cast<std::chrono::microseconds> (loop_end - loop_start).count()) * pow(10, -6);

            std::printf("done: f=%d/%d  in %fs \n", counter, this->n_f_vals, loop_dur);

            f = sqrt(f);
            counter++;
            mc_sweep = 0;
            for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->NM; ++i)
                this->hist[i] = 0;
        }
    }

    this->normalize_JDOS();
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
        for (int i = 0; i < this->ising_lattice->NE; ++i)
            if (ln_JDOS[i * this->ising_lattice->NM + q] > 0)
                temp += exp(this->ln_JDOS[i * this->ising_lattice->NM + q] - this->ln_JDOS[first_idx * this->ising_lattice->NM + q]);

        long double sum_ln_JDOS = ln_JDOS[first_idx * this->ising_lattice->NM + q] + log(temp);

        for (int i = 0; i < this->ising_lattice->NE; ++i)
            if (this->ln_JDOS[i * this->ising_lattice->NM + q] > 0)
                this->ising_lattice->JDOS[i * this->ising_lattice->NM + q] = exp(this->ln_JDOS[i * this->ising_lattice->NM + q] + this->ising_lattice->norm_factor[q] - sum_ln_JDOS);
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


#include "wl.h"

WL::WL(double f_init, double f_final, double flatness, u_int64_t seed) {
    this->init(f_init, f_final, flatness, seed);
}

WL::WL(double f_init, double f_final, double flatness, u_int64_t seed, System &ising_lattice) {
    std::printf("add in the constructor: %p \n", &(ising_lattice));
    this->init(f_init, f_final, flatness, seed);
    this->add_lattice(ising_lattice);
    std::printf("add in the constructor after: %p \n", &(ising_lattice));
}

WL::~WL() {
    std::printf("here wl\n");
    delete this->ising_lattice;
    delete[] this->ln_JDOS;
    delete[] this->hist;
    delete[] this->f_vals;
}

void WL::init(double f_init, double f_final, double flatness, uint64_t seed) {
    this->f_init = f_init;
    this->f_final = f_final;
    this->flatness = flatness;

    this->seed = seed;

    this->n_f_vals = 0;
    double f = f_init;
    while(f > f_final) {
        this->n_f_vals++;
        f = sqrt(f);
    }
    
    this->f_vals = new double[this->n_f_vals];
    f = f_init;
    int i = 0;
    while (f > f_final) {
        this->f_vals[i] = f;
        i++;
        f = sqrt(f);
    }
}

void WL::add_lattice(System &ising_lattice) {
    std::printf("add in the add_latice: %p \n", &(ising_lattice));
    this->ising_lattice = &(ising_lattice);
    std::printf("add in the add_lattice after: %p \n", &(ising_lattice));
    std::printf("add in the add_lattice after (this): %p \n", this->ising_lattice);

    this->ln_JDOS = new long double[this->ising_lattice->NE * this->ising_lattice->NM];
    this->hist = new long long[this->ising_lattice->NE * this->ising_lattice->NM];

    this->ising_lattice->dim=6;
}



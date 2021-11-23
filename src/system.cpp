
#include "system.h"

System::System(int L, int Sz, std::string lattice) {
    this->L = L;
    this->Sz = Sz;

    this->lattice = lattice;

    if (this->lattice == "SS") {
        this->dim = 2;
        this->N_atm = L*L;
        this->NN = 4;
    } else if (this->lattice == "SC") {
        this->dim = 3;
        this->N_atm = L*L*L;
        this->NN = 6;
    } else if (this->lattice == "BCC") {
        this->dim = 3;
        this->N_atm = 2 * L*L*L;
        this->NN = 8;
    } else if (this->lattice == "FCC") {
        this->dim = 3;
        this->N_atm = 4 * L*L*L;
        this->NN = 12;
    } else if (this->lattice == "HCP") {
        this->dim = 3;
        this->N_atm = 2 * L*L*L;
        this-> NN = 12;
    } else if (this->lattice == "HEX") {
        this->dim = 3;
        this->N_atm = L*L*L;
        this->NN = 8;
    } else {
        std::printf(" -- Error invalid Lattice -- \n");
    }

    int max_E = (1.0 / 2.0) * NN * N_atm;
    int max_M = N_atm;

    this->NE = 1 + (max_E / 2.0);
    this->NM = max_M + 1;

    this->energies = this->create_map(- max_E, max_E, 4);
    this->magnetizations = this->create_map(- max_M, max_M, 2);

    std::string NN_table_file_name = "../neighbour_tables/neighbour_table_" + std::to_string(this->dim) + "D_" + this->lattice +
    "_" + std::to_string(this->NN) + "NN_L" + std::to_string(this->L) + ".txt";
    std::string norm_factor_file_name = "../coefficients/coefficients_" + std::to_string(this->N_atm) + "d" + std::to_string(this->Sz) + ".txt";

    this->norm_factor = new long double[this->NM];
    this->NN_table = new int[this->N_atm * this->NN];

    this->read_norm_factor(norm_factor_file_name);
    this->read_NN_talbe(NN_table_file_name);

    this->spins_vector = new int[this->N_atm];
    this->E_config = this->energy();
    this->M_config = this->magnetization();

    this->JDOS = new long double[this->NE * this->NM];
    for (int i = 0; i < this->NE * this->NM; ++i)
        this->JDOS[i] = 0;
}

System::System() {}

System::~System() {
    delete[] this->norm_factor;
    delete[] this->NN_table;
    delete[] this->spins_vector;
    delete[] this->JDOS;
}

void System::init_spins_max_M() {
    for (int i = 0; i < this->N_atm; ++i) {
        this->spins_vector[i] = 1;
    }

    this->E_config = this->energy();
    this->M_config = this->magnetization();
}

void System::init_spins_random(RNG &rng) {
    for (int i = 0; i < this->N_atm; ++i) {
        if (rng.rand_uniform() < 0.5)
            this->spins_vector[i] = 1;
        else
            this->spins_vector[i] = -1;
    }

    this->E_config = this->energy();
    this->M_config = this->magnetization();
}

int System::energy() {
    int E_config = 0;
    for (int i = 0; i < this->N_atm; ++i) {
        for (int a = 0; a < this->NN; ++a)
            E_config += - this->spins_vector[i] *
            this->spins_vector[this->NN_table[i * this->NN + a]];
    }
    return E_config / 2;
}

int System::magnetization() {
    int M_config = 0;
    for (int i = 0; i < this->N_atm; ++i) {
        M_config += spins_vector[i];
    }
    return M_config;
}

int System::energy_flip(int site) {
    int delta_E = 0;
    for (int a = 0; a < this->NN; ++a) {
        delta_E += - this->spins_vector[site] * this->spins_vector[this->NN_table[site * this->NN + a]];
    }
    return delta_E;
}

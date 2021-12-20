
#include "system.h"

System::System(int L, int Sz, std::string lattice, std::string dir) {
    this->L = L;
    this->Sz = Sz;
    this->S = (this->Sz - 1) / 2.0;

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

    this->max_E = 4 * this->S * this->S * this->NN * this->N_atm / 2.0;
    this->max_M = 2 * this->S * N_atm;

    this->NE = 1 + (this->max_E / 2.0);
    this->NM = this->max_M + 1;

    this->energies = this->create_map(- max_E, max_E, 4);
    this->magnetizations = this->create_map(- max_M, max_M, 2);

    std::string NN_table_file_name = dir + "neighbour_tables/neighbour_table_" + std::to_string(this->dim) + "D_" + this->lattice +
    "_" + std::to_string(this->NN) + "NN_L" + std::to_string(this->L) + ".txt";
    std::string norm_factor_file_name = dir + "coefficients/coefficients_" + std::to_string(this->N_atm) + "d" + std::to_string(this->Sz) + ".txt";
    std::string Npos_file_name = dir + "sum_npos/sum_configs_Npos" + std::to_string(this->Sz) + "_N_atm" + std::to_string(this->N_atm) + ".txt";

    this->norm_factor = new long double[this->NM];
    this->NN_table = new int[this->N_atm * this->NN];

    this->read_norm_factor(norm_factor_file_name);
    this->read_NN_talbe(NN_table_file_name);
    this->read_N_pos(Npos_file_name);

    this->spins_vector = new int[this->N_atm];
    
    this->spins_values = new int[this->Sz];
    for (int i = 0; i < this->Sz; ++i) {
        this->spins_values[i] = 2 * this->S - 2 * i;
    }


    this->JDOS = new long double[this->NE * this->NM];
    for (int i = 0; i < this->NE * this->NM; ++i) {
        this->JDOS[i] = 0;
    }
}

System::~System() {
    delete[] this->norm_factor;
    delete[] this->NN_table;
    delete[] this->spins_vector;
    delete[] this->JDOS;
    for (int i = 0; i < this->line_size_Npos.size(); i++)
        delete[] this->Npos[i];
    delete[] this->Npos;
}

void System::init_spins_max_M() {
    for (int i = 0; i < this->N_atm; ++i) {
        this->spins_vector[i] = this->spins_values[0];
    }

    this->E_config = this->energy();
    this->M_config = this->magnetization();
}

void System::init_spins_random(RNG &rng) {
    for (int i = 0; i < this->N_atm; ++i) {
        int idx = rand() % this->Sz;
        this->spins_vector[i] = this->spins_values[idx];
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

int System::energy_flip(int site, int new_spin_idx) {
    int E_tmp1 = 0;
    for (int a = 0; a < this->NN; ++a) {
        E_tmp1 += - this->spins_vector[site] * this->spins_vector[this->NN_table[site * this->NN + a]];
    }
    int E_tmp2 = 0;
    for (int a = 0; a < this->NN; ++a) {
        E_tmp2 += - this->spins_values[new_spin_idx] * this->spins_vector[this->NN_table[site * this->NN + a]];
    }
    return E_tmp2 - E_tmp1;
}

int System::magnetization_flip(int site, int new_spin_idx) {
    return this->spins_values[new_spin_idx] - this->spins_vector[site];
}

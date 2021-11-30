
#ifndef FSS_H
#define FSS_H

#include <cmath>
#include <vector>
#include <array>
#include <stdint.h>

#include "system.h"
#include "rng.h"

class FSS {
    private:
        long long REP;
        int skip;
        int idx_M_max;
        int n_f_vals;

        unsigned long long *hist;
        std::vector<int> *flip_list;
        
        System *ising_lattice;
        RNG *rng;

        double *time_iter;
        unsigned long long *steps_iter;

        bool added_rng = false;
        bool added_params = false;
        bool added_lattice = false;

        unsigned long long min_hist();
        int normalize_JDOS(int &q);
        void first_step();
        void random_config_q(int &q, int &idx_E_config);
        void scan(int &q, int &idx_E_config);
        void spin_flip(int *new_spins_vector, int &new_E_config, int &new_idx_E_config, int &flipped_idx1, int &flipped_idx2);

    public:
        double run_time;

        FSS(RNG &rng);
        FSS(RNG &rng, System &ising_lattice);
        FSS(long long REP, int skip, RNG &rng);
        FSS(long long REP, int skip, RNG &rng, System &ising_lattice);
        ~FSS();

        void set_lattice(System &ising_system);
        void set_rng(RNG &rng);
        void set_params(long long REP, int skip);
        void simulate(int run=0, bool verbose=false);

        void write_to_file(std::string name, bool debug=true);
        void print_JDOS();
};


#endif


#ifndef WL_H
#define WL_H

#include <cmath>
#include <stdint.h>

#include "system.h"
#include "rng.h"

class WL {
    private:
        double f_init;
        double f_final;
        double flatness;
        int n_f_vals;

        long long *hist;
        long double *ln_JDOS;
        long long mc_cycle;
        
        System *ising_lattice;
        RNG *rng;

        double *time_iter;
        long long *steps_iter;

        void init(double, double, double, RNG &);
        long long min_hist();
        long double average_hist();

        void normalize_JDOS();
        bool flat_hist();

    public:
        double run_time;

        WL(double f_init, double f_final, double flatness, RNG &rng);
        WL(double f_init, double f_final, double flatness, RNG &rng, System &ising_lattice);
        ~WL();

        void add_lattice(System &ising_system);
        void simulate(long long steps, int run=0, bool verbose=false);
        void write_to_file(std::string name, bool debug=true);
        
        

};

#endif

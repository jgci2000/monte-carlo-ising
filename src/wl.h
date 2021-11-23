
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

        void init(double, double, double, RNG &);
        long long min_hist();
        long double average_hist();

        void normalize_JDOS();
        bool flat_hist();

    public:
        WL(double, double, double, RNG &);
        WL(double, double, double, RNG &, System &);
        ~WL();

        void add_lattice(System &);
        void simulate(long long);
        
        

};

#endif

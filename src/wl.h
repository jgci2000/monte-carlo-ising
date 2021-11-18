
#ifndef WL_H
#define WL_H

#include <cmath>
#include <stdint.h>

#include "system.h"

class WL {
    private:
        double f_init;
        double f_final;
        double flatness;

        uint64_t seed;

        
        long long *hist;
        long double *ln_JDOS;
        long long mc_cycle;

        double *f_vals;
        int n_f_vals;

        void init(double, double, double, uint64_t);

    public:

        System *ising_lattice;

        WL(double, double, double, uint64_t);
        WL(double, double, double, uint64_t, System &);
        ~WL();

        void add_lattice(System &);

};

#endif

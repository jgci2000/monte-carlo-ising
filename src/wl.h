/**
 * Implementation by João Inácio (j.g.c.inacio@fys.uio.no).
 * Dec. 2021
 */

#ifndef WL_H
#define WL_H

#include <cmath>
#include <stdint.h>

#include "system.h"
#include "rng.h"

/**
 * Wang Landau class. 
 * Implements the Wang Landau algorithm for the calculation 
 * of the Joint Density of States for the Ising model. 
 * A System and RNG objects need to be added to this class so the 
 * algorithm works as intended.
 * Can output the number of steps for each iteration and the time taken.
 * One can access the total run time trough the variable run_time.
 */
class WL {
    private:
        double f_init;
        double f_final;
        double flatness;
        int n_f_vals;

        unsigned long long *hist;
        long double *ln_JDOS;
        unsigned long long mc_cycle;
        
        System *ising_lattice;
        RNG *rng;

        double *time_iter;
        unsigned long long *steps_iter;

        bool added_rng = false;
        bool added_params = false;
        bool added_lattice = false;

        unsigned long long min_hist();
        long double average_hist();

        void normalize_JDOS();
        bool flat_hist();

    public:
        double run_time;

        WL(RNG &rng);
        WL(RNG &rng, System &ising_lattice);
        WL(double f_init, double f_final, double flatness, RNG &rng);
        WL(double f_init, double f_final, double flatness, RNG &rng, System &ising_lattice);
        ~WL();

        void set_lattice(System &ising_system);
        void set_rng(RNG &rng);
        void set_params(double f_init, double f_final, double flatness);
        void simulate(unsigned long long steps, int run=0, bool verbose=false);

        void write_to_file(std::string name, std::string dir, bool debug=true);
        void print_JDOS();
};

#endif

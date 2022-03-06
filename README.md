# Monte-Carlo Methods for the Ising model

Here you can find all of the codes and results used in the paper titled *Accurate Estimate of the Joint Density of States via Flat Scan Sampling* with the reference arXiv:submit/4195427 [cond-mat.stat-mech] 5 Mar 2022.

Any questions contact: <j.g.c.inacio@fys.uio.no>

## Suported Models:
 * Ising Model
 * Ising SpinS Model

## Current Methods:
 * Wang-Landau 
 * Flat Scan Sampling 

## Prerequisits
 * C/C++ compiler (compiler has to support the c++17 standard library). 

# How to use

To use the developed code you have to include `system.h` and `rng.h`. Next, depending on the algorithm you want to use, include either `fss.h` or `wl.h`. The `RNG` class takes as input the seed and uses the xoshiro256** algorithm for a fast generation of random numbers. The `System` class encapsules all of the system's properties and results. It uses as inputs the number unit cells `L`, the spin degeneracy in the z-direction, the lattice type and finally the relative path to where the normalization, neighbor tables and sum_pos files are stored.

For the `FSS` method, it can be used like this:
```cpp
#include <iostream>

#include "../src/system.h"
#include "../src/fss.h"
#include "../src/rng.h"

int main(int argv, char **argc) {
    RNG rng(20);
    System ising_lattice(4, 2, "SS", "../");

    FSS fss(1000, 1, rng, ising_lattice);

    fss.simulate(0, true);
    fss.write_to_file("JDOS.txt", "./", true);
    
    return 0;
}
```
The `FSS` class takes as parameters the parameter REP and skip, an instance of the `RNG` and `System` classes. 

For the `WL` method, it can be used like this:
```cpp
#include <iostream>

#include "../src/system.h"
#include "../src/wl.h"
#include "../src/rng.h"

int main(int argv, char **argc) {
    RNG rng(20);
    System ising_lattice(4, 2, "SS");
    
    WL wang_landau(exp(1), 1 + pow(10, -8), 0.9, rng, ising_lattice);
    
    wang_landau.simulate(10000, 0, true);
    wang_landau.print_JDOS();
    
    return 0;
}
```
The `WL` class takes as parameters the value of the first and last modification factors, the flatness criteria and an instance of `RNG` and `System`. 

For spinS computations, it it better to use the MATLAB implementation by J. Amaral <jamaral@ua.pt>. A more efficient C++ implementation will come in the near future.
		


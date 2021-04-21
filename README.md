# Monte-Carlo Methods in Stat. Mech.

Any questions contact: ```inacio.joao16@ua.pt```

## Current Methods:
 * Wang-Landau (C++)
 * Flat Scan Sampling (C++ and MPI)

## Prerequisits
 * C/C++ compiler, preferably the GNU compiler. 
 * MPI libraries and compiler. 

# Instructions

First off to change the system size and lattice, you ahve to open the source file and chanhe the ´´´#define´´´.

Run the command ```make``` in the respective folder to compile the source code. To compile the MPI version of FSS, use ```make fss_mpi```. This will generate an executable file in the respective directory. 

To run the algorithm just type ```./name_of_executable args``` in the console. To run the MPI version use ```mpirun -np number_cores ./name_of_executable args``` . If no CLI arguments are passed, the computations will be done with default parameters. To change the parameters, you have to pass arguments in the following order:

## WL:
  * 1º -> exponent to f_final (```f_final = 1 + pow(10, - agrument)```). The default value is 8.
  * 2º -> flatness parameter. The default is 90. For large simulations 90-95 is a good flatness criteria. 
  * 3º -> number of the run. Used to numerate multiple runs of the simulation in oder to obtain statistics from the method. Default is 0. 
	
There are 3 cases for providing arguments to the WL sampling simulation. Just the exponent of f_final, the exponent of f_final and run, respectively, and exponent, run and flatness criteria, in this order.

## FSS single core:
  * 1º -> exponent of REP (```REP = pow(10, argument)```). The default is 4.
  * 2º -> skip parameter. The default is N_atm. This number is probably fine, unless you would like to compare simulations with different skip values. 
  * 3º -> number of the run. Default if 0.
	
There are also 3 cases. Just the REP exponent, the REP exponent and run, and the REP exponent, skip and run, in this order.

## FSS MPI:
  * 1º -> exponent of REP. Default is 4.
  * 2º -> exponent shuffle REP, it is computed in the same way as the REP vaule. This number dictates the number of times the spins configuration is shuffled before starting the computation for a given magnetization. Default is the value of REP. A good value may be 100.REP to 1000.REP. This a very quick step, does not affect the simulation time in any way.
  * 3º -> skip. Default if N_atm.
  * 4ª number of the run. Default is 0.

There are 4 cases for giving arguments to the simulation. Just the REP exponent, the REP exponent and run, the REP exponent, shuffle exponent and run, and REP exponent, shuffle exponent, skip and run, in this order. 

## Results:

Once the algorithm is done, you will find three files in the directory ´´´/data´´´ in each of the source folders of both methods. One is the JDOS computed by the method, the other a copy of the prints to the terminal throught the simulation and the last is a file which contains some usefull information about wall times in step of the simulation and total wall time. 




		


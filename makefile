#######################
# make all
#
all: wl_Ising fss_Ising fss_spinS

#######################
# make each one
#
wl_Ising: functions.o
	g++ ./src/wl_Ising.cpp -o wl_Ising -Ofast -std=c++17 -m64 ./include/functions.o

fss_Ising: functions.o
	g++ ./src/fss_Ising.cpp -o fss_Ising -Ofast -std=c++17 -m64 ./include/functions.o
	mpic++ ./src/fss_MPI_Ising.cpp -o fss_MPI_Ising -Ofast -std=c++17 -m64 ./include/functions.o

fss_spinS: functions.o
	g++ ./src/fss_Ising_SpinS.cpp -o fss_spinS -Ofast -std=c++17 -m64 ./include/functions.o

########################
# Dependencies
#
functions.o: ./include/functions.cpp
	g++ -c ./include/functions.cpp -o ./include/functions.o -Ofast -std=c++17 -m64

########################
# make clean
#
clean: clean_wl_Ising clean_fss_Ising clean_fss_spinS clean_functions

clean_wl_Ising: 
	rm wl_Ising

clean_fss_Ising:
	rm fss_Ising
	rm fss_MPI_Ising

clean_fss_spinS:
	rm fss_spinS

clean_functions:
	rm ./include/*.o

clean_results:
	rm -r ./data/


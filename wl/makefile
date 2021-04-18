wl: WL_Functions.o
	g++ wl_Ising.cpp -o wl_Ising -Ofast -std=c++17 -m64 -march=native WL_Functions.o

WL_Functions.o: WL_Functions.cpp
	g++ -c WL_Functions.cpp -o WL_Functions.o

clean_wl:
	rm *.o
	rm wl_Ising

clean_results:
	rm -r ./data/**

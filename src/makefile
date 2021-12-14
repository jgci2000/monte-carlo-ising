
test_wl: test_wl.cpp system.cpp system.h wl.h wl.cpp rng.h 
	g++-11 test_wl.cpp system.cpp wl.cpp -o test_wl

test_fss: test_fss.cpp system.cpp system.h fss.h fss.cpp rng.h 
	g++-11 test_fss.cpp system.cpp fss.cpp -o test_fss

test_spins: test_spinS.cpp system.cpp system.h wl.h wl.cpp rng.h fss.h fss.cpp fss_spinS.cpp 
	g++-11 test_spinS.cpp system.cpp wl.cpp fss.cpp fss_spinS.cpp -o test_spins

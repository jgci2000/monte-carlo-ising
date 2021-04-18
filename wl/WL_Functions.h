//
// Header file for FSS method functions.
// João Inácio, Mar. 25th, 2021
//


#ifndef WL_FUNCTIONS_H
#define WL_FUNCTIONS_H

struct system_info 
{
    int L;
    int N_atm;
    int NN;
    std::string lattice;
    int dim;
};

long long min_hist(long long *, int);

long double average_hist(long long *, int);

std::map<int, int> create_map(int, int, int);

void read_NN_talbe(std::string, int *);

void read_norm_factor(std::string, long double *);

system_info get_system(int, int);

std::vector<std::string> split(const std::string&, char);

#endif // FSS_FUNCTIONS_H

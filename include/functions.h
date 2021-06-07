//
// Header file for FSS method functions.
// João Inácio, Apr. 17th, 2021
//


#ifndef FSS_FUNCTIONS_H
#define FSS_FUNCTIONS_H

struct system_info 
{
    int L;
    int N_atm;
    int NN;
    std::string lattice;
    int dim;
    int SZ;
};

struct xorshift64s_state {
  uint64_t a;
};

long long min_hist(long long *, int);

long double average_hist(long long *, int);

void shuffle(long double *, long long, std::array<std::vector<int>, 2> &, int *, 
int, int, int, int, int *, int &, int &, std::map<int, int> &);

std::map<int, int> create_map(int, int, int);

void read_NN_talbe(std::string, int *);

void read_norm_factor(std::string, long double *, int);

system_info get_system(int, int);

system_info get_system(int, int, double);

std::vector<std::string> split(const std::string&, char);

std::vector<int> split_int(const std::string&, char);

#endif // FSS_FUNCTIONS_H

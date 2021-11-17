//
// Source file for FSS method functions.
// João Inácio, Apr. 17th, 2021
//


#include <iostream>
#include <fstream>
#include <climits>
#include <cmath>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <stdint.h>

#include "functions.h"


uint64_t xorshift64s(struct xorshift64s_state *state)
{
	uint64_t x = state->a;
	x ^= x >> 12;
	x ^= x << 25;
	x ^= x >> 27;
	state->a = x;
	return x * UINT64_C(0x2545F4914F6CDD1D);
}


long long min_hist(long long *hist, int size) 
{
    long long min = __LONG_LONG_MAX__;
    for (int i = 0; i < size; i++) 
        if (hist[i] != 0 && hist[i] < min)
            min = hist[i];
    return min;
}

long double average_hist(long long *hist, int size)
{
    long double sum = 0;
    int nnz = 0;
    for (int i = 0; i < size; i++)
    {
        if (hist[i] != 0)
        {
            sum += hist[i];
            nnz++;
        }
    }
    return sum / nnz;
}

void shuffle(long double *JDOS, long long REP, std::array<std::vector<int>, 2> &flip_list, int *spins_vector, int q, int N_atm, int NN, int NM, int *NN_table, int &E_config, int &idx_E_config, std::map<int, int> &energies)
{
    srand((unsigned) time(NULL));
    xorshift64s_state state = {.a = (uint64_t) rand()};

    int *new_spins_vector = new int[N_atm];
    bool accepted = false;

    // Shuffle for REP times
    for (int idx_shuffle = 0; idx_shuffle < REP; idx_shuffle++)
    {
        // Get a new random condig at magnetization q
        if (!accepted)
            for (int i = 0; i < N_atm; i++)
                new_spins_vector[i] =  spins_vector[i];
        int new_E_config = 0;
        int new_idx_E_config = 0;

        // Flip a positive spin to a negative
        int idx_tmp1 = rand() % flip_list[0].size();
        int flipped_idx1 = flip_list[0].at(idx_tmp1);
        new_spins_vector[flipped_idx1] = - 1;

        flip_list[1].push_back(flipped_idx1);
        flip_list[0].erase(flip_list[0].begin() + idx_tmp1);

        int delta_E = 0;
        for (int a = 0; a < NN; a++)
            delta_E += - new_spins_vector[flipped_idx1] * new_spins_vector[NN_table[flipped_idx1 * NN + a]];
        new_E_config = E_config + 2 * delta_E;

        // Flip a negative spin to a positive
        int idx_tmp2 = rand() % flip_list[1].size();
        int flipped_idx2 = flip_list[1].at(idx_tmp2);
        new_spins_vector[flipped_idx2] = 1;

        flip_list[0].push_back(flipped_idx2);
        flip_list[1].erase(flip_list[1].begin() + idx_tmp2);

        delta_E = 0;
        for (int a = 0; a < NN; a++)
            delta_E += - new_spins_vector[flipped_idx2] * new_spins_vector[NN_table[flipped_idx2 * NN + a]];
        new_E_config = new_E_config + 2 * delta_E;          

        // Wang Landau criteria
        new_idx_E_config = energies[new_E_config];
        long double ratio = JDOS[idx_E_config * NM + q] / JDOS[new_idx_E_config * NM + q];

        if (ratio >= 1 || ((long double) xorshift64s(&state) / (long double) __DBL_MAX__) < ratio)
        {
            for (int i = 0; i < N_atm; i++)
                spins_vector[i] = new_spins_vector[i];

            E_config = new_E_config;
            idx_E_config = new_idx_E_config;
            accepted = true;
        }
        else
        {
            if (flipped_idx1 != flipped_idx2)
            {
                flip_list[0].pop_back();
                flip_list[0].push_back(flipped_idx1);

                flip_list[1].pop_back();
                flip_list[1].push_back(flipped_idx2);
            }
            accepted = false;
        }
    }
}

std::map<int, int> create_map(int init, int final, int step)
{
    std::map<int, int> out;
    int i = 0;
    while (init <= final)
    {
        out.insert(std::pair<int, int>(init, i));
        init += step;
        i++;
    }
    return out;
}

void read_NN_talbe(std::string file_name, int *NN_table)
{
    std::ifstream neighbour_tables_file(file_name);
    std::string line;

    if (neighbour_tables_file.is_open())
    {
        int i = 0;
        while (std::getline(neighbour_tables_file, line))
        {
            std::vector<std::string> a = split(line, ' ');
            for (int idx = 0; idx < a.size(); idx++)
                NN_table[i++] = std::stold(a.at(idx));
        }
        neighbour_tables_file.close();
    }
    else
        std::cout << "Unable to open neighbour table file. Invalid lattice size or lattice type." << std::endl;
}

void read_norm_factor(std::string file_name, long double *norm_factor, int N_atm)
{
    std::ifstream norm_factor_file(file_name);
    std::string line;

    if (N_atm <= 1024)
    {
        if (norm_factor_file.is_open()) 
        {
            for (int i = 0; std::getline(norm_factor_file, line); i++)
                norm_factor[i] = log(std::stold(line));
            norm_factor_file.close();
        }
        else 
            std::cout << "Unable to open normalization factor file. Invalid lattice size or the file isn't on the correct directory." << std::endl;
    }
    else
    {
        if (norm_factor_file.is_open()) 
        {
            for (int i = 0; std::getline(norm_factor_file, line); i++)
                norm_factor[i] = std::stold(line);
            norm_factor_file.close();
        }
        else 
            std::cout << "Unable to open normalization factor file. Invalid lattice size or the file isn't on the correct directory." << std::endl;
    }

}

system_info get_system(int L, int lattice_num)
{
    system_info system;
    system.L = L;

    switch (lattice_num)
    {
        case 1:
            system.dim = 2;
            system.lattice = "SS";
            system.N_atm = L * L;
            system.NN = 4;
            break;

        case 2:
            system.dim = 3;
            system.lattice = "SC";
            system.N_atm = L * L * L;
            system.NN = 6;
            break;

        case 3:
            system.dim = 3;
            system.lattice = "BCC";
            system.N_atm = 2 * L * L * L; 
            system.NN = 8;
            break;
        
        case 4:
            system.dim = 3;
            system.lattice = "FCC";
            system.N_atm = 4 * L * L * L;
            system.NN = 12;
            break;

        case 5: 
            system.dim = 3;
            system.lattice = "HCP";
            system.N_atm = 2 * L * L * L;
            system.NN = 12;
            break;
        case 6:
            system.dim = 3;
            system.lattice = "Hex";
            system.N_atm = L * L * L;
            system.NN = 8;
            break;

        default:
            std::cout << "Invalid lattice number." << std::endl;
            break;
    }

    return system;
}

system_info get_system(int L, int lattice_num, double S)
{
    system_info system;
    system.L = L;
    system.SZ = 2*S + 1;

    switch (lattice_num)
    {
        case 1:
            system.dim = 2;
            system.lattice = "SS";
            system.N_atm = L * L;
            system.NN = 4;
            break;

        case 2:
            system.dim = 3;
            system.lattice = "SC";
            system.N_atm = L * L * L;
            system.NN = 6;
            break;

        case 3:
            system.dim = 3;
            system.lattice = "BCC";
            system.N_atm = 2 * L * L * L; 
            system.NN = 8;
            break;
        
        case 4:
            system.dim = 3;
            system.lattice = "FCC";
            system.N_atm = 4 * L * L * L;
            system.NN = 12;
            break;

        case 5: 
            system.dim = 3;
            system.lattice = "HCP";
            system.N_atm = 2 * L * L * L;
            system.NN = 12;
            break;
        case 6:
            system.dim = 3;
            system.lattice = "Hex";
            system.N_atm = L * L * L;
            system.NN = 8;
            break;

        default:
            std::cout << "Invalid lattice number." << std::endl;
            break;
    }

    return system;
}

std::vector<std::string> split(const std::string& s, char seperator)
{
    std::vector<std::string> output;

    std::string::size_type prev_pos = 0, pos = 0;

    while((pos = s.find(seperator, pos)) != std::string::npos)
    {
        std::string substring( s.substr(prev_pos, pos-prev_pos) );

        output.push_back(substring);

        prev_pos = ++pos;
    }

    output.push_back(s.substr(prev_pos, pos-prev_pos)); // Last word

    return output;
}

std::vector<int> split_int(const std::string& s, char seperator)
{
    std::vector<int> output;

    std::string::size_type prev_pos = 0, pos = 0;

    while((pos = s.find(seperator, pos)) != std::string::npos)
    {
        std::string substring( s.substr(prev_pos, pos-prev_pos) );

        output.push_back(stoi(substring));

        prev_pos = ++pos;
    }

    output.push_back(stoi(s.substr(prev_pos, pos-prev_pos))); // Last word

    return output;
}

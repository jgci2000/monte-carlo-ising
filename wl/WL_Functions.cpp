//
// Source file for FSS method functions.
// João Inácio, Mar. 25th, 2021
//


#include <iostream>
#include <fstream>
#include <climits>
#include <array>
#include <vector>
#include <map>
#include <string>

#include "WL_Functions.h"


long long min_hist(long long *hist, int size) 
{
    long long min = LONG_LONG_MAX;
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

void read_norm_factor(std::string file_name, long double *norm_factor)
{
    std::ifstream norm_factor_file(file_name);
    std::string line;

    if (norm_factor_file.is_open()) 
    {
        for (int i = 0; std::getline(norm_factor_file, line); i++)
            norm_factor[i] = std::stold(line);
        norm_factor_file.close();
    }
    else 
        std::cout << "Unable to open normalization factor file. Invalid lattice size or the file isn't on the correct directory." << std::endl;
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

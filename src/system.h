
#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdint.h>

#include "rng.h"

class System {
    public:
        int L;
        int N_atm;
        int Sz;
        int NN;

        int dim;
        std::string lattice;

        int NE;
        int NM;
        std::map<int, int> energies;
        std::map<int, int> magnetizations;

        int *NN_table;
        long double *norm_factor;

        int *spins_vector;
        long double *JDOS;

        int E_config;
        int M_config;

        System();
        System(int, int, std::string);
        ~System();

        void init_spins_random(RNG &rng);
        void init_spins_max_M();

        int energy();
        int magnetization();

        int energy_flip(int);

    private:
        std::map<int, int> create_map(int init, int final, int step) {
            std::map<int, int> out;
            int i = 0;
            while (init <= final) {
                out.insert(std::pair<int, int>(init, i));
                init += step;
                i++;
            }
            return out;
        }

        std::vector<std::string> split(const std::string& s, char seperator) {
            std::vector<std::string> output;
            std::string::size_type prev_pos = 0, pos = 0;
            while((pos = s.find(seperator, pos)) != std::string::npos) {
                std::string substring( s.substr(prev_pos, pos-prev_pos) );
                output.push_back(substring);
                prev_pos = ++pos;
            }
            output.push_back(s.substr(prev_pos, pos-prev_pos)); // Last word
            return output;
        }

        void read_NN_talbe(std::string file_name) {
            std::ifstream neighbour_tables_file(file_name);
            std::string line;
            if (neighbour_tables_file.is_open()) {
                int i = 0;
                while (std::getline(neighbour_tables_file, line)) {
                    std::vector<std::string> a = split(line, ' ');
                    for (int idx = 0; idx < a.size(); idx++)
                        this->NN_table[i++] = std::stold(a.at(idx));
                }
                neighbour_tables_file.close();
            }
            else
                std::cout << "Unable to open neighbour table file. Invalid lattice size or lattice type." << std::endl;
        }

        void read_norm_factor(std::string file_name) {
            std::ifstream norm_factor_file(file_name);
            std::string line;
            if (this->N_atm <= 1024) {
                if (norm_factor_file.is_open()) {
                    for (int i = 0; std::getline(norm_factor_file, line); i++)
                        this->norm_factor[i] = log(std::stold(line));
                    norm_factor_file.close();
                }
                else
                    std::cout << "Unable to open normalization factor file. Invalid lattice size or the file isn't on the correct directory." << std::endl;
            }
            else {
                if (norm_factor_file.is_open()) {
                    for (int i = 0; std::getline(norm_factor_file, line); i++)
                        this->norm_factor[i] = std::stold(line);
                    norm_factor_file.close();
                }
                else
                    std::cout << "Unable to open normalization factor file. Invalid lattice size or the file isn't on the correct directory." << std::endl;
            }
        }
};

#endif


#include "fss.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <stdint.h>

void FSS::simulate_Sz_S(int run, bool verbose) {
    
    int *new_spins_vector = new int[this->ising_lattice->N_atm];

    std::vector<int> SPM;
    for (int i = 0; i < this->ising_lattice->N_atm; i++)
        SPM.push_back(0);

    int prev_idx_E_config = 0;
    int prev_idx_Npos = 0;

    auto runtime_start = std::chrono::steady_clock::now();

    this->first_step_Sz_S(SPM);

    if (verbose) {
        time_t now = time(0);
        std::string t = ctime(&now); t.pop_back();
        now = time(0);
        t = ctime(&now); t.pop_back();
        printf("%s | run: %d | q: %d/%d \n", 
            t.c_str(),
            run, 
            0, 
            this->idx_M_max);
    }
    
    for (int q = 1; q <= this->idx_M_max; q++)
    {
        auto q_start = std::chrono::steady_clock::now();

        long long *hist2 = new long long[this->ising_lattice->NE * this->ising_lattice->line_size_Npos.at(q)];

        for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->line_size_Npos.at(q); i++)
        {
            hist2[i] = 0;
        }

        // Random config at q

        for (int i = 0; i < this->ising_lattice->N_atm; i++)
            this->ising_lattice->spins_vector[i] = this->ising_lattice->spins_values[0];
        int E_config =  - this->ising_lattice->max_E;
        int idx_E_config = this->ising_lattice->energies[E_config];

        for (int i = 0; i < this->ising_lattice->N_atm; i++)
            SPM[i] = 0;

        for (int idx = 1; idx <= q; idx++)
        {
            std::vector<int> flip_list;
            for (int i = 0; i < this->ising_lattice->N_atm; i++)
                if (SPM[i] + 1 < this->ising_lattice->Sz)
                    flip_list.push_back(i);
            int flipped_idx = flip_list.at(this->rng->rand() % flip_list.size());

            int E_old = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                E_old += - this->ising_lattice->spins_vector[flipped_idx] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flipped_idx * this->ising_lattice->NN + a]];

            SPM[flipped_idx]++;
            this->ising_lattice->spins_vector[flipped_idx] = this->ising_lattice->spins_values[SPM[flipped_idx]];

            int E_new = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                E_new += - this->ising_lattice->spins_vector[flipped_idx] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flipped_idx * this->ising_lattice->NN + a]];

            E_config = E_config - E_old + E_new;
        }

        std::vector<int> counter; counter.assign(this->ising_lattice->Sz, 0);
        for (int i = 0; i < this->ising_lattice->N_atm; i++)
            for (int j = 0; j < this->ising_lattice->Sz; j++)
                if (SPM[i] == j)
                    counter[j]++;

        int counter2 = 0;
        int idx_Npos;
        for (idx_Npos = 0; idx_Npos < this->ising_lattice->line_size_Npos.at(q) * this->ising_lattice->Sz; idx_Npos += this->ising_lattice->Sz)
        {
            for (int i = 0; i < this->ising_lattice->Sz; i++)
                if (counter[i] == this->ising_lattice->Npos[q][i + idx_Npos])
                    counter2++;

            if (counter2 == this->ising_lattice->Sz)
                break;
            else
                counter2 = 0;
        }

        idx_E_config = this->ising_lattice->energies[E_config];
        int idx_Npos_config = idx_Npos / this->ising_lattice->Sz;

        // Update hist22ograms

        hist2[idx_E_config * this->ising_lattice->line_size_Npos.at(q) + idx_Npos_config]++;

        // Scan the first config
        // printf("HERE_before! \n");

        for (int x = 1; x < this->ising_lattice->Sz; x++)
        {
            std::vector<int> flip_list;
            for (int i = 0; i < this->ising_lattice->N_atm; i++)
                if (SPM[i] + 1 <= this->ising_lattice->Sz - x)
                    flip_list.push_back(i);

            // for (int i = 0; i < flip_list.size(); ++i)
            //     printf("%d ", flip_list.at(i));
            // printf("\n");

            for (int flip_idx = 0; flip_idx < flip_list.size(); flip_idx++)
            {
                std::vector<int> SPM_tmp = SPM;

                int E_tmp1 = 0;
                for (int a = 0; a < this->ising_lattice->NN; a++)
                    E_tmp1 += - this->ising_lattice->spins_vector[flip_list.at(flip_idx)] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flip_list.at(flip_idx) * this->ising_lattice->NN + a]];

                int E_tmp2 = 0;
                for (int a = 0; a < this->ising_lattice->NN; a++)
                    E_tmp2 += - this->ising_lattice->spins_values[SPM[flip_list.at(flip_idx)] + x] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flip_list.at(flip_idx) * this->ising_lattice->NN + a]];

                int E_tmp3 = E_config - E_tmp1 + E_tmp2;
                int idx_E_tmp3 = this->ising_lattice->energies[E_tmp3];

                SPM_tmp[flip_list.at(flip_idx)] = SPM[flip_list.at(flip_idx)] + x;

                std::vector<int> counter; counter.assign(this->ising_lattice->Sz, 0);
                for (int i = 0; i < this->ising_lattice->N_atm; i++)
                    for (int j = 0; j < this->ising_lattice->Sz; j++)
                        if (SPM_tmp[i] == j)
                            counter[j]++;

                int counter2 = 0;
                int idx_Npos;
                for (idx_Npos = 0; idx_Npos < this->ising_lattice->line_size_Npos.at(q + x) * this->ising_lattice->Sz; idx_Npos += this->ising_lattice->Sz)
                {
                    for (int i = 0; i < this->ising_lattice->Sz; i++)
                        if (counter[i] == this->ising_lattice->Npos[q + x][i + idx_Npos])
                            counter2++;

                    if (counter2 == this->ising_lattice->Sz)
                        break;
                    else
                        counter2 = 0;
                }
                idx_Npos /= this->ising_lattice->Sz;

                this->JDOS_M_spin[q + x][idx_Npos * this->ising_lattice->NE + idx_E_tmp3] += this->JDOS_M_spin[q][idx_Npos_config * this->ising_lattice->NE + idx_E_config] / this->REP;
            }
        }

        // for (int i = 0; i < 3; ++i) {
        //     for (int j = 0; j < this->ising_lattice->line_size_Npos.at(i) * this->ising_lattice->NE; ++j) {
        //         printf("%f ", this->JDOS_M_spin[i][j]);
        //     }
        //     printf("\n");
        // }


        long long k = 1;
        bool accepted = false;
        std::vector<int> new_SPM;

        // Where the magic happens

        while (min_hist2(hist2, this->ising_lattice->NE * this->ising_lattice->line_size_Npos.at(q)) < this->REP)
        {
            // Get a new random condig at magnetization q
            // printf("HERE! \n");
            if (!accepted)
            {
                new_SPM = SPM;
                for (int i = 0; i < this->ising_lattice->N_atm; i++)
                    new_spins_vector[i] =  this->ising_lattice->spins_vector[i];
            }
            int new_E_config = 0;
            int new_idx_E_config = 0;
            int new_idx_Npos_config = 0;

            // Choose a spin to flip

            int flipped_idx1 = this->rng->rand() % this->ising_lattice->N_atm;

            int E_old = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                E_old += - new_spins_vector[flipped_idx1] * new_spins_vector[this->ising_lattice->NN_table[flipped_idx1 * this->ising_lattice->NN + a]];

            std::vector<int> SPM_end_list;
            for (int i = 0; i < this->ising_lattice->Sz; i++)
                if (i != SPM[flipped_idx1])
                    SPM_end_list.push_back(i);

            int SPM_start = SPM[flipped_idx1];
            int SPM_end = SPM_end_list.at(this->rng->rand() % SPM_end_list.size());

            new_spins_vector[flipped_idx1] = this->ising_lattice->spins_values[SPM_end];
            new_SPM[flipped_idx1] = SPM_end;

            int E_new = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                E_new += - new_spins_vector[flipped_idx1] * new_spins_vector[this->ising_lattice->NN_table[flipped_idx1 * this->ising_lattice->NN + a]];
            new_E_config = E_config - E_old + E_new;

            int SPM_dif = SPM_end - SPM_start;

            // Choose another one to flip back

            std::vector<int> flip_list;
            for (int i = 0; i < this->ising_lattice->N_atm; i++)
                if (new_SPM[i] + 1 - SPM_dif >= 1 && new_SPM[i] + 1 - SPM_dif <= this->ising_lattice->Sz)
                    flip_list.push_back(i);

            int flipped_idx2 = flip_list.at(this->rng->rand() % flip_list.size());

            E_old = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                E_old += - new_spins_vector[flipped_idx2] * new_spins_vector[this->ising_lattice->NN_table[flipped_idx2 * this->ising_lattice->NN + a]];

            new_SPM[flipped_idx2] = new_SPM[flipped_idx2] - SPM_dif;
            new_spins_vector[flipped_idx2] = this->ising_lattice->spins_values[new_SPM[flipped_idx2]];

            E_new = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                E_new += - new_spins_vector[flipped_idx2] * new_spins_vector[this->ising_lattice->NN_table[flipped_idx2 * this->ising_lattice->NN + a]];
            new_E_config = new_E_config - E_old + E_new;

            std::vector<int> counter; counter.assign(this->ising_lattice->Sz, 0);
            for (int i = 0; i < this->ising_lattice->N_atm; i++)
                for (int j = 0; j < this->ising_lattice->Sz; j++)
                    if (new_SPM[i] == j)
                        counter[j]++;

            int counter2 = 0;
            int idx_Npos;
            for (idx_Npos = 0; idx_Npos < this->ising_lattice->line_size_Npos.at(q) * this->ising_lattice->Sz; idx_Npos += this->ising_lattice->Sz)
            {
                for (int i = 0; i < this->ising_lattice->Sz; i++)
                    if (counter[i] == this->ising_lattice->Npos[q][i + idx_Npos])
                        counter2++;

                if (counter2 == this->ising_lattice->Sz)
                    break;
                else
                    counter2 = 0;
            }

            new_idx_Npos_config = idx_Npos / this->ising_lattice->Sz;
            new_idx_E_config = this->ising_lattice->energies[new_E_config];

            // Wang Landau criteria

            long double ratio = this->JDOS_M_spin[q][idx_Npos_config * this->ising_lattice->NE + idx_E_config] / this->JDOS_M_spin[q][new_idx_Npos_config * this->ising_lattice->NE + new_idx_E_config];

            if (ratio >= 1 || this->rng->rand_uniform() < ratio || hist2[new_idx_E_config * this->ising_lattice->line_size_Npos.at(q) + new_idx_Npos_config] == 0)
            {
                for (int i = 0; i < this->ising_lattice->N_atm; i++)
                    this->ising_lattice->spins_vector[i] = new_spins_vector[i];

                E_config = new_E_config;
                idx_E_config = new_idx_E_config;
                idx_Npos_config = new_idx_Npos_config;
                SPM = new_SPM;

                accepted = true;
            }
            else
            {
                accepted = false;
            }
            // printf("HERE 2! \n");

            // Scan configuration

            if (hist2[idx_E_config * this->ising_lattice->line_size_Npos.at(q) + idx_Npos_config] < this->REP && k % this->skip == 0 || hist2[idx_E_config * this->ising_lattice->line_size_Npos.at(q) + idx_Npos_config] == 0)
            {
                for (int x = 1; x < this->ising_lattice->Sz; x++)
                {
                    std::vector<int> flip_list;
                    for (int i = 0; i < this->ising_lattice->N_atm; i++)
                        if (SPM[i] + 1 <= this->ising_lattice->Sz - x)
                            flip_list.push_back(i);

                    for (int flip_idx = 0; flip_idx < flip_list.size(); flip_idx++)
                    {
                        std::vector<int> SPM_tmp = SPM;

                        int E_tmp1 = 0;
                        for (int a = 0; a < this->ising_lattice->NN; a++)
                            E_tmp1 += - this->ising_lattice->spins_vector[flip_list.at(flip_idx)] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flip_list.at(flip_idx) * this->ising_lattice->NN + a]];

                        int E_tmp2 = 0;
                        for (int a = 0; a < this->ising_lattice->NN; a++)
                            E_tmp2 += - this->ising_lattice->spins_values[SPM[flip_list.at(flip_idx)] + x] * this->ising_lattice->spins_vector[this->ising_lattice->NN_table[flip_list.at(flip_idx) * this->ising_lattice->NN + a]];

                        int E_tmp3 = E_config - E_tmp1 + E_tmp2;
                        int idx_E_tmp3 = this->ising_lattice->energies[E_tmp3];

                        SPM_tmp[flip_list.at(flip_idx)] = SPM[flip_list.at(flip_idx)] + x;

                        std::vector<int> counter; counter.assign(this->ising_lattice->Sz, 0);
                        for (int i = 0; i < this->ising_lattice->N_atm; i++)
                            for (int j = 0; j < this->ising_lattice->Sz; j++)
                                if (SPM_tmp[i] == j)
                                    counter[j]++;

                        int counter2 = 0;
                        int idx_Npos;
                        for (idx_Npos = 0; idx_Npos < this->ising_lattice->line_size_Npos.at(q + x) * this->ising_lattice->Sz; idx_Npos += this->ising_lattice->Sz)
                        {
                            for (int i = 0; i < this->ising_lattice->Sz; i++)
                                if (counter[i] == this->ising_lattice->Npos[q + x][i + idx_Npos])
                                    counter2++;

                            if (counter2 == this->ising_lattice->Sz)
                                break;
                            else
                                counter2 = 0;
                        }
                        idx_Npos /= this->ising_lattice->Sz;

                        this->JDOS_M_spin[q + x][idx_Npos * this->ising_lattice->NE + idx_E_tmp3] += this->JDOS_M_spin[q][idx_Npos_config * this->ising_lattice->NE + idx_E_config] / this->REP;
                    }
                }

                hist2[idx_E_config * this->ising_lattice->line_size_Npos.at(q) + idx_Npos_config]++;
            }
        
            k++;
            // printf("%d \n", k);
        }
        
        // Normalize JDOS and output to console

        std::vector<long double> sum_JDOS_M_spin; sum_JDOS_M_spin.assign(this->ising_lattice->NE, 0);
        long double sum_sum_JDOS_M_spin = 0;
        for (int i = 0; i < this->ising_lattice->NE; i++)
        {
            for (int j = 0; j < this->ising_lattice->line_size_Npos.at(q + 1); j++)
                sum_JDOS_M_spin[i] += this->JDOS_M_spin[q + 1][j * this->ising_lattice->NE + i];
            sum_sum_JDOS_M_spin += sum_JDOS_M_spin[i];
        }

        for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->line_size_Npos.at(q + 1); i++)
            this->JDOS_M_spin[q + 1][i] = this->JDOS_M_spin[q + 1][i] * exp(this->ising_lattice->norm_factor[q + 1]) / sum_sum_JDOS_M_spin;

        for (int i = 0; i < this->ising_lattice->NE; i++)
            sum_JDOS_M_spin[i] = 0;
        sum_sum_JDOS_M_spin = 0;
        for (int i = 0; i < this->ising_lattice->NE; i++)
        {
            for (int j = 0; j < this->ising_lattice->line_size_Npos.at(q + 1); j++)
                sum_JDOS_M_spin[i] += this->JDOS_M_spin[q + 1][j * this->ising_lattice->NE + i];
            sum_sum_JDOS_M_spin += sum_JDOS_M_spin[i];
        }

        for (int i = 0; i < this->ising_lattice->NE; i++)
            this->ising_lattice->JDOS[i * this->ising_lattice->NM + q + 1] = sum_JDOS_M_spin[i] * exp(this->ising_lattice->norm_factor[q + 1]) / sum_sum_JDOS_M_spin;

        int hits = 0;
        for (int i = 0; i < this->ising_lattice->NE * this->ising_lattice->line_size_Npos.at(q); i++)
            if (this->JDOS_M_spin[q][i] > 0)
                hits++;

        auto q_end = std::chrono::steady_clock::now();
        double q_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (q_end - q_start).count()) * pow(10, -6);

        if (verbose) {
            time_t now = time(0);
            std::string t = ctime(&now); t.pop_back();
            now = time(0);
            t = ctime(&now); t.pop_back();
            printf("%s | run: %d | q: %d/%d | time: %fs | E: %ld | time/E: %fs \n", 
                t.c_str(),
                run, 
                q, 
                this->idx_M_max, 
                q_time, 
                hits, 
                q_time / hits);
        }

        delete[] hist2;
    }

    // {
    //     this->time_iter[q - 1] = q_time;
    //     this->steps_iter[q - 1] = k;
    // }

    auto runtime_end = std::chrono::steady_clock::now();
    this->run_time = (double) (std::chrono::duration_cast<std::chrono::microseconds> (runtime_end - runtime_start).count()) * pow(10, -6);

    delete[] new_spins_vector;
}

void FSS::first_step_Sz_S(std::vector<int> &SPM) {
    this->JDOS_M_spin[0][0] = 1;
    this->ising_lattice->JDOS[0] = 1;

    for (int x = 1; x < this->ising_lattice->Sz; x++)
    {
        for (int flip_idx = 0; flip_idx < this->ising_lattice->N_atm; flip_idx++)
        {
            std::vector<int> SPM_tmp = SPM;

            int E_tmp1 = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                E_tmp1 += - this->ising_lattice->spins_values[SPM[flip_idx]] * this->ising_lattice->spins_values[SPM[flip_idx]];

            int E_tmp2 = 0;
            for (int a = 0; a < this->ising_lattice->NN; a++)
                E_tmp2 += - this->ising_lattice->spins_values[SPM[flip_idx] + x] * this->ising_lattice->spins_values[SPM[flip_idx]];

            int E_tmp3 = - this->ising_lattice->max_E - E_tmp1 + E_tmp2;
            int idx_E_tmp3 = this->ising_lattice->energies[E_tmp3];

            SPM_tmp[flip_idx] = SPM[flip_idx] + x;

            std::vector<int> counter; counter.assign(this->ising_lattice->Sz, 0);
            for (int i = 0; i < this->ising_lattice->N_atm; i++)
                for (int j = 0; j < this->ising_lattice->Sz; j++)
                    if (SPM_tmp[i] == j)
                        counter[j]++;

            int counter2 = 0;
            int idx_Npos;
            for (idx_Npos = 0; idx_Npos < this->ising_lattice->line_size_Npos.at(x) * this->ising_lattice->Sz; idx_Npos += this->ising_lattice->Sz)
            {
                for (int i = 0; i < this->ising_lattice->Sz; i++)
                    if (counter[i] == this->ising_lattice->Npos[x][i + idx_Npos])
                        counter2++;

                if (counter2 == this->ising_lattice->Sz)
                    break;
                else
                    counter2 = 0;
            }
            idx_Npos /= this->ising_lattice->Sz;

            this->JDOS_M_spin[x][idx_Npos * this->ising_lattice->NE + idx_E_tmp3] += this->JDOS_M_spin[0][0];
        }
    }

    std::vector<long double> sum_JDOS_M_spin; sum_JDOS_M_spin.assign(this->ising_lattice->NE, 0);
    long double sum_sum_JDOS_M_spin = 0;
    for (int i = 0; i < this->ising_lattice->NE; i++)
    {
        for (int j = 0; j < this->ising_lattice->line_size_Npos.at(1); j++)
            sum_JDOS_M_spin[i] += this->JDOS_M_spin[1][j * this->ising_lattice->NE + i];
        sum_sum_JDOS_M_spin += sum_JDOS_M_spin[i];
    }

    for (int i = 0; i < this->ising_lattice->NE; i++)
        this->ising_lattice->JDOS[i * this->ising_lattice->NM + 1] = sum_JDOS_M_spin[i] * exp(this->ising_lattice->norm_factor[1]) / sum_sum_JDOS_M_spin;
}



long double average_hist22(long long *hist2, int size)
{
    long double sum = 0;
    int nnz = 0;
    for (int i = 0; i < size; i++)
    {
        if (hist2[i] != 0)
        {
            sum += hist2[i];
            nnz++;
        }
    }
    return sum / nnz;
}

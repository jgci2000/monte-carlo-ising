%
% Explicit calcs for DOS at q=4 for Ising spin=1 2x2 lattice spin 1 (Npos 3)
%
% v1 - works, but fraction of q = 3, E = 0 is a manual input, which can be obtained by scanning, which is OK!
% v2 - FSS for q = 4
% v3 - generalizing for any q and Npos
% v4 - cleanup, QoL (output, temp_JDOS)
% v5 - new k_skip implementation
% v6 - newer skip implementation, skip++ only on accepted step
%
clear
close all
%
L = 8;
Npos = 5; % (2*S)+1
%
dim = '2D';
lattice_type = 'SS';
neighbours = '1NN';
%
N_atm = L^2;
NN = 4;
%
REP = 1E1;
skip = 1E2;
%
JDOS_filename = ['JDOS_FSS_spinS_Ising_v6_Npos', int2str(Npos), '_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip)];
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of ', 'FSS_spinS_Ising_v6_Npos', int2str(Npos), '_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip)]);
%
eval(['load ./coefficients/coefficients_', int2str(N_atm), 'd', int2str(Npos),'.mat'])
eval(['load ./neighbour_tables/neighbour_table_', dim, '_', lattice_type, '_', neighbours, '_L', int2str(L), '.mat'])
%
Z_spin_values(:,1) = -(Npos-1) : 2 : (Npos-1); % always integers
%
M_list(:,1) = -N_atm*max(Z_spin_values) : 2 : N_atm*max(Z_spin_values);
E_list(:,1) = max(Z_spin_values(:,1)).^2*(- N_atm * NN ./2) : 4 : max(Z_spin_values(:,1)).^2*(N_atm * NN ./2); % possible energy values
%
JDOS_nRPS = zeros(length(E_list), length(M_list));
JDOS_nRPS(1,1) = 1;
JDOS_nRPS(1,length(M_list)) = 1;
%
q_max = (length(M_list)-1)/2; %(length(M_list)-2);
%
if (q_max == (length(M_list)-2))
    %
    JDOS_filename = horzcat(JDOS_filename, '_full');
    %
end
%
output = nan(q_max, 5);
%
t_total = tic;
%
for q = 1:q_max
    %
    q_timer = tic;
    %
    hist_nRPS = zeros(length(E_list), 1);
    hist_E_selected_nRPS = zeros(length(E_list), 1);
    %
    % RANDOM SPIN CONFIGURATION AT q
    [S_vector, SPM, E] = ...
        function_random_spin_config_at_q_REDUX(N_atm, NN, Z_spin_values, q, Npos, NN_table);
    %
    hist_nRPS(E_list == E, 1) = hist_nRPS(E_list == E, 1) + 1;
    hist_E_selected_nRPS(E_list == E, 1) = hist_E_selected_nRPS(E_list == E, 1) + 1;
    %
    % SCAN spin S
    %
    for step_type = 1:(Npos-1)
        %
        for i = 1:N_atm
            %
            if SPM(i) <= Npos - step_type
                %
                E_temp1 = -S_vector(i)*( S_vector(NN_table(i,1)) + S_vector(NN_table(i,2)) + S_vector(NN_table(i,3)) + S_vector(NN_table(i,4)) );
                %
                S_vector(i) = Z_spin_values(SPM(i) + step_type);
                %
                E_temp2 = -S_vector(i)*( S_vector(NN_table(i,1)) + S_vector(NN_table(i,2)) + S_vector(NN_table(i,3)) + S_vector(NN_table(i,4)) );
                %
                E_scan = E - E_temp1 + E_temp2;
                JDOS_nRPS(E_list == E_scan, q+step_type) = JDOS_nRPS(E_list == E_scan, q+step_type) + JDOS_nRPS(E_list == E, q);
                %
                S_vector(i) = Z_spin_values(SPM(i));
                %
            end
            %
        end
        %
    end
    %
    k = 1;
    k_skip = 0;
    %
    while min(hist_E_selected_nRPS(hist_E_selected_nRPS > 0)) < REP % CHECK FOR FULL CONFIG SET
        %
        % RANDOM WALK TRIAL STEP AT q
        [S_vector_new, SPM_new, E_new] = ...
            function_rw_step_at_q_REDUX(S_vector, SPM, N_atm, Npos, E, NN_table, Z_spin_values);
        %
        % ACCEPT/REJECT WITH WL WEIGH FACTOR
        %
        if rand < min( [ JDOS_nRPS(E_list == E, q) ./ JDOS_nRPS(E_list == E_new, q), 1]) % accept
            %
            E = E_new;
            S_vector = S_vector_new;
            SPM = SPM_new;
            %
            hist_nRPS(E_list == E_new, 1) = hist_nRPS(E_list == E_new, 1) + 1;
            k_skip = k_skip + 1;
            %
        else % reject
            %
            hist_nRPS(E_list == E, 1) = hist_nRPS(E_list == E, 1) + 1;
            %
        end
        %
        % SCAN
        if (hist_E_selected_nRPS(E_list == E, 1) < REP && k_skip >= skip) || hist_E_selected_nRPS(E_list == E, 1) == 0
            %
            for step_type = 1:(Npos-1)
                %
                for i = 1:N_atm
                    %
                    if SPM(i) <= Npos - step_type
                        %
                        E_temp1 = -S_vector(i)*( S_vector(NN_table(i,1)) + S_vector(NN_table(i,2)) + S_vector(NN_table(i,3)) + S_vector(NN_table(i,4)) );
                        %
                        S_vector(i) = Z_spin_values(SPM(i) + step_type);
                        %
                        E_temp2 = -S_vector(i)*( S_vector(NN_table(i,1)) + S_vector(NN_table(i,2)) + S_vector(NN_table(i,3)) + S_vector(NN_table(i,4)) );
                        %
                        E_scan = E - E_temp1 + E_temp2;
                        JDOS_nRPS(E_list == E_scan, q+step_type) = JDOS_nRPS(E_list == E_scan, q+step_type) + JDOS_nRPS(E_list == E, q);
                        %
                        S_vector(i) = Z_spin_values(SPM(i));
                        %
                    end
                    %
                end
                %
            end
            %
            hist_E_selected_nRPS(E_list == E, 1) = hist_E_selected_nRPS(E_list == E, 1) + 1;
            k_skip = 0;
            %
        end
        %
        k = k + 1;
        %
    end
    %
    JDOS_nRPS(:,q+1) = JDOS_nRPS(:,q+1) ./ sum(JDOS_nRPS(:,q+1)) * norm_factor(q+1);
    %
    hits = nnz(JDOS_nRPS(:,q));
    q_timer = toc(q_timer);
    %
    output(q,:) = [q, q_timer, hits, q_timer/hits, k];
    %
    eval(['save ', JDOS_filename, '_temp_q', int2str(q), '.mat JDOS_nRPS E_list M_list output -v7.3'])
    eval(['delete ', JDOS_filename, '_temp_q', int2str(q-1), '.mat'])
    %
    disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(q), '/', int2str(q_max), ' | time: ', num2str(q_timer), ' secs | E pts: ', int2str(hits), ' | time per E pt: ', num2str(q_timer/hits), ' secs | rw steps: ', int2str(k)])
    %
end
%
JDOS_nRPS(1,length(M_list)) = 1;
%
t_total = toc(t_total);
%
disp([datestr(now,'dd-mm-yyyy HH:MM:SS'), ' finish | total run time ', num2str(t_total), ' seconds'])
%
if (q_max == (length(M_list)-1)/2)
    %
    JDOS_nRPS(:, q_max + 2: q_max + Npos) = zeros(length(E_list), length(1:(Npos-1)));
    %
end
%
eval(['save ', JDOS_filename, '.mat JDOS_nRPS E_list M_list output -v7.3'])
eval(['delete ', JDOS_filename, '_temp_q', int2str(q_max), '.mat'])
%
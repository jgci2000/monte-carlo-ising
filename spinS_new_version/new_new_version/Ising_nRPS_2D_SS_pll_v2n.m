clear
close all
%
% v1 - from serial v1 (which comes from spinS_v6)
% v2 - include hist_E_Mqp1 (failsafe DOS = 0 can't be used in parfor)
%
rng default
%
L = 8;
%
% q_max = 9;
%
dim = '2D';
lattice_type = 'SS';
neighbours = '1NN';
%
REP_total = 1E2; % number of desired configurations per (E,M) pair
skip = L^2; % 1 for no skip
%
n_cores = 10;
% parpool(n_cores);
%
% END OF USER INPUT
%
REP = REP_total/n_cores;
%
N_atm = L^2;
NN = 4;
%
if rem(length(-N_atm : 2 : N_atm),2) == 0
    %
    q_max = (length(-N_atm : 2 : N_atm))/2 - 1;
    %
else
    %
    q_max = (length(-N_atm : 2 : N_atm) + 1)/2 - 1; 
    %
end
%
disp(['Ising ', dim, ' ', lattice_type, ' ', neighbours])
disp(['L = ', int2str(L)])
disp(['REP_total = 1E', num2str(log10(REP_total))])
disp(['skip = ', int2str(skip)])
disp(['n_cores = ', int2str(n_cores)])
%
wspace_filename = ['workspace_nRPS_Ising_pll_v2_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP_total)), '_skip_', int2str(skip), '_ncores_', int2str(n_cores)];
JDOS_filename = ['JDOS_nRPS_Ising_pll_v2_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP_total)), '_skip_', int2str(skip), '_ncores_', int2str(n_cores)];
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | start run of ', 'nRPS_Ising_v72_', dim, '_', lattice_type, '_L', int2str(L), '_REP_1E', int2str(log10(REP)), '_skip_', int2str(skip), '_ncores_', int2str(n_cores)]);
%
eval(['load ./coefficients/coefficients_', int2str(N_atm), 'd', int2str(2),'.mat'])
eval(['load ./neighbour_tables/neighbour_table_', dim, '_', lattice_type, '_', neighbours, '_L', int2str(L), '.mat'])
%
output = nan(q_max, 4); %5);
%
M_list(:,1) = -N_atm : 2 : N_atm;
E_list(:,1) = (- N_atm * NN /2) : 4 : (N_atm * NN /2); % possible energy values
%
JDOS_nRPS = zeros(length(E_list), length(M_list));
JDOS_nRPS(1,1) = 1;
JDOS_nRPS(1,length(M_list)) = 1;
%
t_total = tic;
%
% SCAN AT q = 1, ADD TO JDOS AT q = 2
[JDOS_nRPS] = ...
    function_Ising_scan_norm_correct_q1(N_atm, E_list, NN, NN_table, JDOS_nRPS);
%
% CALC JDOS at q = 2
JDOS_nRPS(:,2) = JDOS_nRPS(:,2) ./ sum(JDOS_nRPS(:,2)) .* norm_factor(2);
%
output(1,:) = [1, 0, 1, 0, ]; %0];
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(1), '/', int2str(q_max)]);
%
% MAIN LOOP
%
for q = 2:q_max
    %
    q_timer = tic;
    %
    JDOS_temp = cell(1, n_cores);
    %
    parfor core = 1:n_cores
        

        %
        JDOS_temp{1, core} = zeros(length(E_list), 1);
        %
        hist_nRPS = zeros(length(E_list), 1);
        hist_E_selected_nRPS = zeros(length(E_list), 1);
        hist_E_Mqp1 = zeros(length(E_list), 1);
        %
        % RANDOM SPIN CONFIGURATION AT q
        [S_vector, E] = ...
            function_Ising_random_spin_config_at_q(N_atm, NN, q, NN_table);
        %
        hist_nRPS(E_list == E, 1) = hist_nRPS(E_list == E, 1) + 1;
        hist_E_selected_nRPS(E_list == E, 1) = hist_E_selected_nRPS(E_list == E, 1) + 1;
        %
        % SCAN
        % [JDOS_nRPS] = ...
        %    function_Ising_scan_norm_correct(S_vector, E, E_list, NN_table, REP, JDOS_nRPS, q);
        %
        for i = 1:N_atm
            %
            if S_vector(i) == 1
                %
                E_temp1 = -S_vector(i)*( S_vector(NN_table(i,1)) + S_vector(NN_table(i,2)) + S_vector(NN_table(i,3)) + S_vector(NN_table(i,4)) );
                %
                E_scan = E - 2*E_temp1;
                hist_E_Mqp1(E_list == E_scan, 1) = hist_E_Mqp1(E_list == E_scan, 1) + 1;
                JDOS_temp{1, core}(E_list == E_scan, 1) = JDOS_temp{1, core}(E_list == E_scan, 1) + JDOS_nRPS(E_list == E, q);
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
            [S_vector_new, E_new] = ...
                function_Ising_rw_step_at_q(S_vector, N_atm, E, NN_table);
            %
            % FAILSAFE FOR DOS = 0
            %
            % if JDOS_nRPS(E_list == E_new, q) == 0
                % %
                % JDOS_nRPS(E_list == E_new, q) = (1/hist_nRPS(E_list == E)) * JDOS_nRPS(E_list == E, q);
                % disp('FAILSAFE FOR DOS = 0')
                % %
            % end            
            %
            % ACCEPT/REJECT WITH WL WEIGH FACTOR
            %
            if rand < min( [ JDOS_nRPS(E_list == E, q) ./ JDOS_nRPS(E_list == E_new, q), 1]) % accept
                %
                E = E_new;
                S_vector = S_vector_new;
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
            %
            if (hist_E_selected_nRPS(E_list == E, 1) < REP && k_skip >= skip) || hist_E_selected_nRPS(E_list == E, 1) == 0
                %
                for i = 1:N_atm
                    %
                    if S_vector(i) == 1
                        %
                        E_temp1 = -S_vector(i)*( S_vector(NN_table(i,1)) + S_vector(NN_table(i,2)) + S_vector(NN_table(i,3)) + S_vector(NN_table(i,4)) );
                        %
                        E_scan = E - 2*E_temp1;
                        hist_E_Mqp1(E_list == E_scan, 1) = hist_E_Mqp1(E_list == E_scan, 1) + 1;
                        JDOS_temp{1, core}(E_list == E_scan, 1) = JDOS_temp{1, core}(E_list == E_scan, 1) + JDOS_nRPS(E_list == E, q);
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
        if min(hist_E_Mqp1(hist_E_Mqp1 > 0)) < REP
            %
            disp('warning : hist_E_Mqp1 below REP - more scanning!')
            %
            REP_temp = REP;
            %
            while min(hist_E_Mqp1(hist_E_Mqp1 > 0)) < REP
                %
                REP_temp = REP_temp + REP;
                %
                while min(hist_E_selected_nRPS(hist_E_selected_nRPS > 0)) < REP_temp % CHECK FOR FULL CONFIG SET
                    %
                    % RANDOM WALK TRIAL STEP AT q
                    [S_vector_new, E_new] = ...
                        function_Ising_rw_step_at_q(S_vector, N_atm, E, NN_table);
                    %
                    % FAILSAFE FOR DOS = 0
                    %
                    % if JDOS_nRPS(E_list == E_new, q) == 0
                        % %
                        % JDOS_nRPS(E_list == E_new, q) = (1/hist_nRPS(E_list == E)) * JDOS_nRPS(E_list == E, q);
                        % disp('FAILSAFE FOR DOS = 0')
                        % %
                    % end
                    %
                    % ACCEPT/REJECT WITH WL WEIGH FACTOR
                    %
                    if rand < min( [ JDOS_nRPS(E_list == E, q) ./ JDOS_nRPS(E_list == E_new, q), 1]) % accept
                        %
                        E = E_new;
                        S_vector = S_vector_new;
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
                    if (hist_E_selected_nRPS(E_list == E, 1) < REP_temp && k_skip >= skip) || hist_E_selected_nRPS(E_list == E, 1) == 0
                        %
                        for i = 1:N_atm
                            %
                            if S_vector(i) == 1
                                %
                                E_temp1 = -S_vector(i)*( S_vector(NN_table(i,1)) + S_vector(NN_table(i,2)) + S_vector(NN_table(i,3)) + S_vector(NN_table(i,4)) );
                                %
                                E_scan = E - 2*E_temp1;
                                hist_E_Mqp1(E_list == E_scan, 1) = hist_E_Mqp1(E_list == E_scan, 1) + 1;
                                JDOS_temp{1, core}(E_list == E_scan, 1) = JDOS_temp{1, core}(E_list == E_scan, 1) + JDOS_nRPS(E_list == E, q);
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
            end
            %
            disp(['min(hist_E_Mqp1) now >= REP, REP_temp = ', int2str(REP_temp)])
            %
        end

   

    %
    end

    for core = 1:n_cores
        %
        JDOS_nRPS(:,q+1) = JDOS_nRPS(:,q+1) + JDOS_temp{1, core}(:, 1);
        %
    end

    %
    JDOS_nRPS(:,q+1) = JDOS_nRPS(:,q+1) ./ sum(JDOS_nRPS(:,q+1)) .* norm_factor(q+1);
    %
    hits = nnz(JDOS_nRPS(:,q));
    q_timer = toc(q_timer);
    %
    output(q,:) = [q, q_timer, hits, q_timer/hits]; %, k];
    %
    eval(['save ', JDOS_filename, '_temp_q', int2str(q), '.mat JDOS_nRPS E_list M_list output -v7.3'])
    eval(['delete ', JDOS_filename, '_temp_q', int2str(q-1), '.mat'])
    %
    disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(q), '/', int2str(q_max), ' | time: ', num2str(q_timer), ' secs | E pts: ', int2str(hits), ' | time per E pt: ', num2str(q_timer/hits), ' secs ']) %| rw steps: ', int2str(k)])
    %
end
%
t_total = toc(t_total);
%
disp([datestr(now,'dd-mm-yyyy HH:MM:SS'), ' finish | total run time ', num2str(t_total), ' seconds'])
% eval(['save ', wspace_filename, '.mat -v7.3'])
eval(['save ', JDOS_filename, '.mat JDOS_nRPS E_list M_list output -v7.3'])
eval(['delete ', JDOS_filename, '_temp_q', int2str(q_max), '.mat'])
%
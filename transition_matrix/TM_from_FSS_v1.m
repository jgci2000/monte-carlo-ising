clear
close all
%
% load TM2.mat
%
% v1 - from flat_E_v7
%
rng default
%
tic
%
L = 4;
REP = 1E3; % FSS REP
skip = L^2;
%
% END OF USER INPUT
%
% SYSTEM PARAMETERS
%
NN = 4;
N_atm = L^2;
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = -1/2 * N_atm * NN : +4 : 1/2 * N_atm * NN; % possible energy values
%
TM = zeros(length(E_list).^2, length(M_list).^2);
% TM(1,1) = 1;
%
eval(['load coefficients_', int2str(N_atm),'d2.mat'])
%
JDOS_aprox = zeros(length(E_list), length(M_list));
JDOS_aprox(1,1) = 1;
log_JDOS_aprox = log(JDOS_aprox); 
%
% VECTOR NEIGHBOUR TABLE
%
nnxpos=nan(L^2,1);
nnxneg=nan(L^2,1);
nnypos=nan(L^2,1);
nnyneg=nan(L^2,1);
%
for i=1:L
    for j=1:L
        [nnxpos(j+(i-1)*L),nnxneg(j+(i-1)*L),nnypos(j+(i-1)*L),nnyneg(j+(i-1)*L)]=function_NN_list_2D_SS(L,i,j);
    end
end
%
q = 1; % ONLY SCAN, USE ALL SPINS DOWN CONFIGS TO CALCULATE DOS(q=2)
%
neo_previous = zeros(length(E_list), length(E_list)); % old new
%
S_vector = -ones(N_atm+1, 1);
E_old = -1/2 * N_atm * NN;
%
for k = 1:REP
    %
    for flipped_pos = 1:N_atm
        %
        delta_E = - 1 .* ( ...
            S_vector(nnxpos(flipped_pos),1) + ...
            S_vector(nnxneg(flipped_pos),1) + ...
            S_vector(nnypos(flipped_pos),1) + ...
            S_vector(nnyneg(flipped_pos),1)) ; % energy of bonds to NN
        %
        E_new = E_old + 2*delta_E;
        %
        neo_previous(E_list == E_old, E_list == E_new) = neo_previous(E_list == E_old, E_list == E_new) + 1;
        %
    end
    %
end
%
% CALCULATE DOS at q = 2
%
for x = find(neo_previous)
    %
    [u,v] = ind2sub(length(E_list),x);
    JDOS_aprox(v, q+1) = neo_previous(x)/sum(sum(neo_previous));
    %
end
%
JDOS_aprox(:, q+1) = JDOS_aprox(:, q+1) ./sum(JDOS_aprox(:, q+1)) .* norm_factor(q+1);
%
neo_previous_trans = neo_previous';
find_neo_previous = find(neo_previous);
%
for x = 1:length(find_neo_previous)
    %
    [u,v] = ind2sub(length(E_list),find_neo_previous(x));
    %
    TM((q-1)*length(M_list)+u, (q)*length(M_list)+v) = neo_previous(u,v) ./ sum(sum(neo_previous)) ./ JDOS_aprox(v, q+1);
    TM((q)*length(M_list)+v, (q-1)*length(M_list)+u) = neo_previous_trans(v,u) ./ sum(sum(neo_previous)) ./ JDOS_aprox(u, q);
    %
end
%
disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: 1/', num2str((length(M_list)-1)/2 + 1 - 1)])
%
for q = 2:N_atm%((length(M_list)-1)/2 + 1 - 1)
    %
    q_time = tic;
    %
    neo_previous = zeros(length(E_list), length(E_list));
    %
    DOS = JDOS_aprox(:, q);
    log_DOS = log(JDOS_aprox(:, q));
    %
    % INITITAL CONFIGURATION
    %
    length_up = (q-1); % N_atm/2;
    % length_down = N_atm - length_up; % N_atm/2;
    %
    site_list(:,1) = randperm(N_atm);
    S_vector = ones(N_atm,1);
    S_vector(site_list(length_up + 1 : end)) = -1;
    %
    % INITIAL ENERGY CALCULATION
    %
    E_old = 0;
    %
    for i=1:L
        for j=1:L
            E_old = E_old - 1/2*S_vector(j+(i-1)*L,1)*(S_vector(nnxpos(j+(i-1)*L))+S_vector(nnxneg(j+(i-1)*L))+S_vector(nnypos(j+(i-1)*L))+S_vector(nnyneg(j+(i-1)*L)));
        end
    end
    %
    % SAME M RANDOM WALK
    %
    hist_E(:,1) = E_list;
    hist_E(:,2) = zeros(length(E_list),1);
    %
    skip_counter = 0;
    %
    while sum(hist_E(:,2)) < nnz(DOS)*REP
        %
        skip_counter = skip_counter + 1;
        %
        flip_down_index = randi(q-1); % randi(N_atm/2);
        flip_up_index = length_up + randi(N_atm - length_up); % length_up + randi(N_atm/2);
        %
        S_vector(site_list(flip_down_index),1) = -1;
        %
        delta_E = + 1 .* ( ...
            S_vector(nnxpos(site_list(flip_down_index)),1) + ...
            S_vector(nnxneg(site_list(flip_down_index)),1) + ...
            S_vector(nnypos(site_list(flip_down_index)),1) + ...
            S_vector(nnyneg(site_list(flip_down_index)),1)) ; % energy of bonds to NN
        %
        E_new = E_old + 2*delta_E; % build the energy matrix
        %
        S_vector(site_list(flip_up_index),1) = 1;
        %
        delta_E = - 1 .* ( ...
            S_vector(nnxpos(site_list(flip_up_index)),1) + ...
            S_vector(nnxneg(site_list(flip_up_index)),1) + ...
            S_vector(nnypos(site_list(flip_up_index)),1) + ...
            S_vector(nnyneg(site_list(flip_up_index)),1)) ; % energy of bonds to NN
        %
        E_new = E_new + 2*delta_E;
        %
        if rand < min( [ exp(log_DOS(E_list == E_old, 1) - log_DOS(E_list == E_new, 1)), 1]) % accept
            %
            site_list([flip_down_index, flip_up_index]) = site_list([flip_up_index, flip_down_index]);
            E_old = E_new;
            %
            if skip_counter >= skip
                %
                if hist_E(E_list == E_new,2) < REP || hist_E(E_list == E_new,2) == 0
                    %
                    hist_E(E_list == E_new,2) = hist_E(E_list == E_new,2) + 1;
                    %
                    % SCAN FOR WL ACCEPT
                    %
                    pos_scan = site_list((length_up+1):end,1);
                    %
                    for spin_index = 1:length(pos_scan)
                        %
                        flipped_pos_scan = pos_scan(spin_index, 1);
                        %
                        delta_E = - 1 .* ( ...
                            S_vector(nnxpos(flipped_pos_scan),1) + ...
                            S_vector(nnxneg(flipped_pos_scan),1) + ...
                            S_vector(nnypos(flipped_pos_scan),1) + ...
                            S_vector(nnyneg(flipped_pos_scan),1)) ; % energy of bonds to NN
                        %
                        neo_previous(E_list == E_old, E_list == (E_old + 2*delta_E)) = neo_previous(E_list == E_old, E_list == (E_old + 2*delta_E)) + 1;
                        %
                    end
                    %
                    skip_counter = 0;
                    %
                end
                %
            end
            %
        else % reject
            %
            S_vector(site_list(flip_down_index),1) = 1;
            S_vector(site_list(flip_up_index),1) = -1;
            %
            if skip_counter >= skip
                %
                if hist_E(E_list == E_old,2) < REP || hist_E(E_list == E_new,2) == 0
                    %
                    hist_E(E_list == E_old,2) = hist_E(E_list == E_old,2) + 1;
                    %
                    % SCAN FOR WL REJECT
                    %
                    pos_scan = site_list((length_up+1):end,1);
                    %
                    for spin_index = 1:length(pos_scan)
                        %
                        flipped_pos_scan = pos_scan(spin_index, 1);
                        %
                        delta_E = - 1 .* ( ...
                            S_vector(nnxpos(flipped_pos_scan),1) + ...
                            S_vector(nnxneg(flipped_pos_scan),1) + ...
                            S_vector(nnypos(flipped_pos_scan),1) + ...
                            S_vector(nnyneg(flipped_pos_scan),1)) ; % energy of bonds to NN
                        %
                        neo_previous(E_list == E_old, E_list == (E_old + 2*delta_E)) = neo_previous(E_list == E_old, E_list == (E_old + 2*delta_E)) + 1;
                        %
                    end
                    %
                    skip_counter = 0;
                    %
                end
                %
            end
            %
        end
        %
    end
    %
    % CALCULATE DOS at q = q+1
    %
    rehits = find(neo_previous);
    %
    for x = 1:length(rehits)
        %
        [u,v] = ind2sub(length(E_list),rehits(x));
        JDOS_aprox(v, q+1) = JDOS_aprox(v, q+1) + neo_previous(rehits(x))./sum(neo_previous(u,:)).*DOS(u,1);
        %
    end
    %
    JDOS_aprox(:, q+1) = JDOS_aprox(:, q+1) ./sum(JDOS_aprox(:, q+1)) .* norm_factor(q+1);
    %
    %
    neo_previous_trans = neo_previous';
    find_neo_previous = find(neo_previous);
    %
    for x = 1:length(find_neo_previous)
        %
        [u,v] = ind2sub(length(E_list),find_neo_previous(x));
        %
        TM((q-1)*length(M_list)+u, (q)*length(M_list)+v) = neo_previous(u,v) ./ REP ./ N_atm / (JDOS_aprox(v, q+1)/JDOS_aprox(u, q));
        TM((q)*length(M_list)+v, (q-1)*length(M_list)+u) = neo_previous_trans(v,u) ./ REP ./ N_atm ;
        %
    end
    %
    q_time = toc(q_time);
    %
    disp([datestr(now,'dd/mm/yyyy HH:MM:SS'), ' | q: ', int2str(q),'/', num2str((length(M_list)-1)/2 + 1 - 1), ' | q_time: ', num2str(q_time), ' | E: ', num2str(nnz(JDOS_aprox(:, q))), ' | q_time/E: ', num2str(q_time/nnz(JDOS_aprox(:, q)))])
    %
end
%
[V,D] = eig(TM);
[~,ind] = sort(diag(D),'descend');
%
JDOS_from_TM = reshape(V(:,ind(1)), [length(E_list), length(M_list)]);
%
for q = 1:(N_atm+1)
    %
    JDOS_from_TM(:,q) = JDOS_from_TM(:,q) ./sum(JDOS_from_TM(:,q)) .* norm_factor(q);
    %
end
%
toc

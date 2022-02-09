clear
clc
close

% Script to visualize the JDoS of an Ising system
% 17th Jan 2022

rep_exp = 3;
L = 4;
S = 1/2;
Sz = uint8(2 * S + 1);
N_atm = L * L;
lattice = "SS";
NN = 4;

max_E = 4 * S * S * NN * N_atm / 2;
max_M = 2 * S * N_atm;

E_list = - max_E:4:max_E;
M_list = - max_M:2:max_M;

file_name = "JDOS_L" + int2str(L) + "_" + lattice + "_Sz_" + int2str(Sz) + "_R1E" + int2str(rep_exp);
JDOS = importdata("spinS/data/" + file_name + ".txt");
index_M0 = (length(M_list)-1)/2 + 1;
JDOS(:,index_M0+1:length(M_list)) = JDOS(:,index_M0-1:-1:1);
z_label = "JDOS";

JDOS = log(JDOS);
z_label = "log(JDOS)";

% code for top view
figure(1)

b = bar3(JDOS);
zlabel('JDOS')
view(0, 90) % top
%view(133,33)
axis([0 length(M_list)+1 0 length(E_list)+1])
xticks([1 ((length(M_list)-1)/4 + 1) ((length(M_list)-1)/2+1) (3*(length(M_list)-1)/4+1) length(M_list)])
yticks([1 ((length(E_list)-1)/4 + 1) ((length(E_list)-1)/2+1) (3*(length(E_list)-1)/4+1) length(E_list)])
%
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
%
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'XTickLabel', linspace(M_list(min(xt)), M_list(max(xt)), numel(xt)))
xlabel('M')%, 'FontSize', 16)
set(gca, 'YTickLabel', linspace(E_list(min(yt)), E_list(max(yt)), numel(yt)))
ylabel('E')%, 'FontSize', 16)
%
h = get(gca,'children');
%
for w = 1:length(h)
    %
    hc = get(h(w),'cdata');
    hz = get(h(w),'zdata');
    %
    for u = 1:(length(hc(:,1))/6)
        %
        if sum(nansum( hc(1+(u-1)*6 : u*6, 1:4) )) == 0
            %
            hc_new = hc;
            hc_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'cdata', hc_new);
            hc = hc_new;
            %
            hz_new = hz;
            hz_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'zdata', hz_new);
            hz = hz_new;
            %
        end
        %
    end
    %
end
%
cbh = colorbar;
cbh.Label.String = z_label;
%cbh.FontSize = 16;
set(gca,'Fontsize',16)

saveas(gcf, "spinS/" + file_name, 'epsc')
rep_exp = 3;
L = 4;
S = 1/2;
Sz = uint8(2 * S + 1);
N_atm = L * L;
lattice = "SS";
NN = 4;

max_E = 4 * S * S * NN * N_atm / 2;
max_M = 2 * S * N_atm;

E_list = - max_E:4:max_E;
M_list = - max_M:2:max_M;

file_name = "JDOS_L" + int2str(L) + "_" + lattice + "_Sz_" + int2str(Sz) + "_R1E" + int2str(rep_exp);
JDOS = importdata("spinS/data/" + file_name + ".txt");
index_M0 = (length(M_list)-1)/2 + 1;
JDOS(:,index_M0+1:length(M_list)) = JDOS(:,index_M0-1:-1:1);
z_label = "JDOS";

JDOS = log(JDOS);
z_label = "log(JDOS)";

% code for top view
figure(1)

b = bar3(JDOS);
zlabel('JDOS')
view(0, 90) % top
%view(133,33)
axis([0 length(M_list)+1 0 length(E_list)+1])
xticks([1 ((length(M_list)-1)/4 + 1) ((length(M_list)-1)/2+1) (3*(length(M_list)-1)/4+1) length(M_list)])
yticks([1 ((length(E_list)-1)/4 + 1) ((length(E_list)-1)/2+1) (3*(length(E_list)-1)/4+1) length(E_list)])
%
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
%
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'XTickLabel', linspace(M_list(min(xt)), M_list(max(xt)), numel(xt)))
xlabel('M')%, 'FontSize', 16)
set(gca, 'YTickLabel', linspace(E_list(min(yt)), E_list(max(yt)), numel(yt)))
ylabel('E')%, 'FontSize', 16)
%
h = get(gca,'children');
%
for w = 1:length(h)
    %
    hc = get(h(w),'cdata');
    hz = get(h(w),'zdata');
    %
    for u = 1:(length(hc(:,1))/6)
        %
        if sum(nansum( hc(1+(u-1)*6 : u*6, 1:4) )) == 0
            %
            hc_new = hc;
            hc_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'cdata', hc_new);
            hc = hc_new;
            %
            hz_new = hz;
            hz_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'zdata', hz_new);
            hz = hz_new;
            %
        end
        %
    end
    %
end
%
cbh = colorbar;
cbh.Label.String = z_label;
%cbh.FontSize = 16;
set(gca,'Fontsize',16)

saveas(gcf, "spinS/" + file_name, 'epsc')

% saveas(gcf,'JDOS_exact_L4_SS.jpeg')


% code for side view
figure(2)

b = bar3(JDOS);
zlabel('JDOS')
% view(0, 90) % top
view(133,33)
axis([0 length(M_list)+1 0 length(E_list)+1])
xticks([1 ((length(M_list)-1)/4 + 1) ((length(M_list)-1)/2+1) (3*(length(M_list)-1)/4+1) length(M_list)])
yticks([1 ((length(E_list)-1)/4 + 1) ((length(E_list)-1)/2+1) (3*(length(E_list)-1)/4+1) length(E_list)])
%
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
%
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'XTickLabel', linspace(M_list(min(xt)), M_list(max(xt)), numel(xt)))
xlabel('M')%, 'FontSize', 16)
set(gca, 'YTickLabel', linspace(E_list(min(yt)), E_list(max(yt)), numel(yt)))
ylabel('E')%, 'FontSize', 16)
%
h = get(gca,'children');
%
for w = 1:length(h)
    %
    hc = get(h(w),'cdata');
    hz = get(h(w),'zdata');
    %
    for u = 1:(length(hc(:,1))/6)
        %
        if sum(nansum( hc(1+(u-1)*6 : u*6, 1:4) )) == 0
            %
            hc_new = hc;
            hc_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'cdata', hc_new);
            hc = hc_new;
            %
            hz_new = hz;
            hz_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'zdata', hz_new);
            hz = hz_new;
            %
        end
        %
    end
    %
end
%
cbh = colorbar;
cbh.Label.String = z_label;
%cbh.FontSize = 16;
set(gca,'Fontsize',16)


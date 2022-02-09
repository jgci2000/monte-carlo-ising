clear
clc
close

% Script to visualize the JDoS of an Ising system
% 17th Jan 2022

save_name = "JDOS_L4_SS";
plot_idx = 0;
title_plot = ["(a) S=1/2", "(b) S=1", "(c) S=3/2", "(d) S=2"];

for S = [1/2, 1, 3/2, 2]
    plot_idx = plot_idx + 1;
    
    rep_exp = 4;
    L = 4;
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
    
    max_JDOS = max(max(JDOS));
    JDOS = JDOS / max_JDOS;

    % code for top view
    figure(1)
    set(gcf,'position',[0, 0, 10000, 4000])
    %fig = gcf;
    %fig.PaperUnits = 'centimeters';
    %fig.PaperPosition = [0 0 12 18];
    subplot(1, 4, plot_idx);

    b = bar3(JDOS);
    zlabel('JDOS')
    view(0, 90)
    %axis([-1 1 -1 1])
    axis square
    xlim([0 max_M + 2])
    ylim([0 max_E * 0.575])
    
    %set(gca, 'Position', [-0.1,0.1,1.2,0.85])
    %set(gcf, 'OuterPosition', [0,0,1,1])
    
    title(title_plot(plot_idx))
       
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
    set(gca, 'XTickLabel', linspace(M_list(min(xt))/max_M, M_list(max(xt))/max_M, numel(xt)))
    xlabel('M')
    set(gca, 'YTickLabel', linspace(E_list(min(yt))/max_E, E_list(max(yt))/max_E, numel(yt)))
    ylabel('E')
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
    set(gca,'Fontsize',15)
end

exportgraphics(gcf, "spinS/" + save_name + ".eps", 'Resolution', 200)


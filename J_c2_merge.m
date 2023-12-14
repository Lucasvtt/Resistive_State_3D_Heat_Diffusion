%--------------------------------------------------------------------------
% \section{1. Determina as correntes críticas I_{c2}
%             para \kappa = 1, T_0 = 0.84T_c e l_z = 3\xi(0),
%             T_0 = 0.75T_c e l_z = 2.4\xi(0)}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos para T_0 = 0.84T_c}
    %----------------------------------------------------------------------

    load 0.84/T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.15/IV.dat
    
    N = size(IV,1);
    I_c2(1) = IV(N-1,2);
    h_s(1) = 0.15;

    load 0.84/T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.25/IV.dat
    
    N = size(IV,1);
    I_c2(2) = IV(N-1,2);
    h_s(2) = 0.25;
    
    load 0.84/T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.35/IV.dat
    
    N = size(IV,1);
    I_c2(3) = IV(N-1,2);
    h_s(3) = 0.35;

    load 0.84/T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.40/IV.dat
    
    N = size(IV,1);
    I_c2(4) = IV(N-1,2);
    h_s(4) = 0.40;    

    load 0.84/T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.45/IV.dat
    
    N = size(IV,1);
    I_c2(5) = IV(N-1,2);
    h_s(5) = 0.45;    

    load 0.84/T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.50/IV.dat
    
    N = size(IV,1);
    I_c2(6) = IV(N-1,2);
    h_s(6) = 0.50;

    load 0.84/T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.75/IV.dat
    
    N = size(IV,1);
    I_c2(7) = IV(N-1,2);
    h_s(7) = 0.75;

    load 0.84/T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_1.00/IV.dat
    
    N = size(IV,1);
    I_c2(8) = IV(N-1,2);
    h_s(8) = 1.00;

    [p,s] = polyfit(h_s(1:8),I_c2(1:8),4);
    h_s_val = 0.15:0.01:1;
    I_c2_val = polyval(p,h_s_val);

    R2 = 1 - (s.normr/norm(I_c2(1:8) - mean(I_c2(1:8))))^2;

    h_s_0p84 = h_s;
    h_s_val_0p84 = h_s_val;
    I_c2_val_0p84 = I_c2_val;
    I_c2_0p84 = I_c2;
    
    clear h_s;
    clear h_val;
    clear I_c2;
    clear I_c2_val;

    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos para T_0 = 0.75T_c}
    %----------------------------------------------------------------------

    load 0.75/T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.15/IV.dat
    
    N = size(IV,1);
    I_c2(1) = IV(N-1,2);
    h_s(1) = 0.15;    

    load 0.75/T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.25/IV.dat
    
    N = size(IV,1);
    I_c2(2) = IV(N-1,2);
    h_s(2) = 0.25;

    load 0.75/T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.30/IV.dat
    
    N = size(IV,1);
    I_c2(3) = IV(N-1,2);
    h_s(3) = 0.30;
    
    load 0.75/T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.50/IV.dat
    
    N = size(IV,1);
    I_c2(4) = IV(N-1,2);
    h_s(4) = 0.50;

    load 0.75/T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.75/IV.dat
    
    N = size(IV,1);
    I_c2(5) = IV(N-1,2);
    h_s(5) = 0.75;

    load 0.75/T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/IV.dat
    
    N = size(IV,1);
    I_c2(6) = IV(N-1,2);
    h_s(6) = 1.00;

    [p,s] = polyfit(h_s(1:6),I_c2(1:6),4);
    h_s_val = 0.15:0.01:1;
    I_c2_val = polyval(p,h_s_val);

    R2 = 1 - (s.normr/norm(I_c2(1:6) - mean(I_c2(1:6))))^2;

    h_s_0p75 = h_s;
    h_s_val_0p75 = h_s_val;
    I_c2_val_0p75 = I_c2_val;
    I_c2_0p75 = I_c2;    

    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico da correte critica em funcao de h_s}
    %----------------------------------------------------------------------
    
    figure;
    l_0p75 = plot(h_s_val_0p75,I_c2_val_0p75,'Color',[0 0.4470 0.7410],'LineWidth',1.5);
    hold on;
    plot(h_s_0p75,I_c2_0p75,'o','Color',[0 0.4470 0.7410],'LineWidth',1.5);
    l_0p84 = plot(h_s_val_0p84,I_c2_val_0p84,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
    plot(h_s_0p84,I_c2_0p84,'o','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
    xlim([min([min(h_s_val_0p75) min(h_s_val_0p84)]) max([max(h_s_val_0p75) max(h_s_val_0p84)])]);
    ylim([0.95*min([min(I_c2_val_0p75) min(I_c2_val_0p84)]) max([max(I_c2_val_0p75) max(I_c2_val_0p84)])]);
    box on;
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$h_s$','Interpreter','LaTeX','FontSize',14);
    ylabel('$I_{c2}/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    pos = 'best';
    leg = legend([l_0p75,l_0p84],'$T_0=0.75T_c,\,l_z=2.4\xi(0)$','$T_0=0.84T_c,\,l_z=3\xi(0)$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('I_c2','_nt.png');
    file_t = strcat('I_c2','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
        
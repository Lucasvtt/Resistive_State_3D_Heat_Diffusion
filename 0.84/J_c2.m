%--------------------------------------------------------------------------
% \section{1. Determina as correntes críticas I_{c2}
%             para \kappa = 1, T_0 = 0.84T_c e l_z = 3\xi(0)}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.15/IV.dat
    
    N = size(IV,1);
    I_c2(1) = IV(N-1,2);
    h_s(1) = 0.15;

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.25/IV.dat
    
    N = size(IV,1);
    I_c2(2) = IV(N-1,2);
    h_s(2) = 0.25;
    
    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.35/IV.dat
    
    N = size(IV,1);
    I_c2(3) = IV(N-1,2);
    h_s(3) = 0.35;

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.40/IV.dat
    
    N = size(IV,1);
    I_c2(4) = IV(N-1,2);
    h_s(4) = 0.40;    

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.45/IV.dat
    
    N = size(IV,1);
    I_c2(5) = IV(N-1,2);
    h_s(5) = 0.45;    

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.50/IV.dat
    
    N = size(IV,1);
    I_c2(6) = IV(N-1,2);
    h_s(6) = 0.50;

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.75/IV.dat
    
    N = size(IV,1);
    I_c2(7) = IV(N-1,2);
    h_s(7) = 0.75;

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_1.00/IV.dat
    
    N = size(IV,1);
    I_c2(8) = IV(N-1,2);
    h_s(8) = 1.00;

    [p,s] = polyfit(h_s(1:8),I_c2(1:8),4);
    h_s_val = 0.15:0.01:1;
    I_c2_val = polyval(p,h_s_val);

    R2 = 1 - (s.normr/norm(I_c2(1:8) - mean(I_c2(1:8))))^2;

    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico da correte critica em funcao de h_s}
    %----------------------------------------------------------------------
    
    figure;
    plot(h_s,I_c2,'o','LineWidth',1.5);
    hold on;
    plot(h_s_val,I_c2_val,'LineWidth',1.5);
    xlim([min(min(h_s)) max(h_s)]);
    ylim([min(min(I_c2)) 1.01*max(I_c2)]);
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$h_s$','Interpreter','LaTeX','FontSize',14);
    ylabel('$I_{c2}/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(0.5,1.565,'$T_0=0.84T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(0.5,1.545,'$l_z=3\xi(0)$','FontSize',24,'Interpreter','LaTeX');
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
    
    cmd = ['mv -f *.png', ' ', 'T_0.84/kappa_1.0/Lz_3'];
    system(cmd);
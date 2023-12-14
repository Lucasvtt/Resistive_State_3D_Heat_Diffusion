%--------------------------------------------------------------------------
% \section{1. Comparacao entre os substratos 
%             para \kappa = 1, T_0 = 0.75T_c e l_z = 2.4\xi(0)}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.25/period.dat
    
    dy = 0.4;
    dz = 0.4;
    l_y = 16;
    l_z = 2.4;
    l_y = l_y-dy;
    l_z = l_z-dz;
    A = l_y*l_z;
    
    J_0p25 = period(:,1); 
    T_0p25 = period(:,2);
    
    J_0p25 = period(:,1);
    T_0p25 = period(:,2);
    
    N = size(J_0p25,1);
    
    JJ_0p25(1) = J_0p25(1);
    m = 1;
    
    for l = 2:N
        if J_0p25(l) == JJ_0p25(m) && J_0p25(l-1) == JJ_0p25(m)
            TT_0p25(m) = T_0p25(l)-T_0p25(l-1);
            F_0p25(m) = 1/TT_0p25(m);
            JJ_0p25(m) = J_0p25(l);
        else
            m = m+1;
            JJ_0p25(m) = J_0p25(l);
        end
    end
    
    I_0p25 = A*JJ_0p25;
    
    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.50/period.dat
    
    J_0p50 = period(:,1);
    T_0p50 = period(:,2);
    
    N = size(J_0p50,1);
    
    JJ_0p50(1) = J_0p50(1);
    m = 1;
    
    for l = 2:N
        if J_0p50(l) == JJ_0p50(m) && J_0p50(l-1) == JJ_0p50(m)
            TT_0p50(m) = T_0p50(l)-T_0p50(l-1);
            F_0p50(m) = 1/TT_0p50(m);
            JJ_0p50(m) = J_0p50(l);
        else
            m = m+1;
            JJ_0p50(m) = J_0p50(l);
        end
    end
    
    I_0p50 = A*JJ_0p50;

    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.75/period.dat
    
    J_0p75 = period(:,1);
    T_0p75 = period(:,2);
    
    N = size(J_0p75,1);
    
    JJ_0p75(1) = J_0p75(1);
    m = 1;
    
    for l = 2:N
        if J_0p75(l) == JJ_0p75(m) && J_0p75(l-1) == JJ_0p75(m)
            TT_0p75(m) = T_0p75(l)-T_0p75(l-1);
            F_0p75(m) = 1/TT_0p75(m);
            JJ_0p75(m) = J_0p75(l);
        else
            m = m+1;
            JJ_0p75(m) = J_0p75(l);
        end
    end
    
    I_0p75 = A*JJ_0p75;

    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/period.dat
    
    J_1p00 = period(:,1);
    T_1p00 = period(:,2);
    
    N = size(J_1p00,1);
    
    JJ_1p00(1) = J_1p00(1);
    m = 1;
    
    for l = 2:N
        if J_1p00(l) == JJ_1p00(m) && J_1p00(l-1) == JJ_1p00(m)
            TT_1p00(m) = T_1p00(l)-T_1p00(l-1);
            F_1p00(m) = 1/TT_1p00(m);
            JJ_1p00(m) = J_1p00(l);
        else
            m = m+1;
            JJ_1p00(m) = J_1p00(l);
        end
    end
    
    I_1p00 = A*JJ_1p00;

    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as frequencias}
    %----------------------------------------------------------------------
    
    figure;
    l_0p25 = plot(I_0p25,F_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,F_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,F_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,F_1p00,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([min([min(F_0p25) min(F_0p50) min(F_0p75) min(F_1p00)]) max([max(F_0p25) max(F_0p50) max(F_0p75) max(F_1p00)])]);
    pos = 'best';
    leg = legend([l_0p25,l_0p50,l_0p75,l_1p00],'$h_s=0.25$','$h_s=0.50$','$h_s=0.75$','$h_s=1.00$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$\nu t_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.05,0.9*max(F_1p00),'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.05,0.83*max(F_1p00),'$l_z=2.4\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('frq','_nt.png');
    file_t = strcat('frq','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0/Lz_2.4'];
    system(cmd);

%--------------------------------------------------------------------------
% \section{2. Comparacao entre os substratos 
%             para \kappa = 1, T_0 = 0.75T_c e l_z = 4.8\xi(0)}
%--------------------------------------------------------------------------

    clear all;
%     close all;
    
    %----------------------------------------------------------------------
    % \subsection{2.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_0.25/period.dat
    
    dy = 0.4;
    dz = 0.4;
    l_y = 16;
    l_z = 4.8;
    l_y = l_y-dy;
    l_z = l_z-dz;
    A = l_y*l_z;
    
    J_0p25 = period(:,1); 
    T_0p25 = period(:,2);
    
    J_0p25 = period(:,1);
    T_0p25 = period(:,2);
    
    N = size(J_0p25,1);
    
    JJ_0p25(1) = J_0p25(1);
    m = 1;
    
    for l = 2:N
        if J_0p25(l) == JJ_0p25(m) && J_0p25(l-1) == JJ_0p25(m)
            TT_0p25(m) = T_0p25(l)-T_0p25(l-1);
            F_0p25(m) = 1/TT_0p25(m);
            JJ_0p25(m) = J_0p25(l);
        else
            m = m+1;
            JJ_0p25(m) = J_0p25(l);
        end
    end
    
    I_0p25 = A*JJ_0p25;
    
    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_0.50/period.dat
    
    J_0p50 = period(:,1);
    T_0p50 = period(:,2);
    
    N = size(J_0p50,1);
    
    JJ_0p50(1) = J_0p50(1);
    m = 1;
    
    for l = 2:N
        if J_0p50(l) == JJ_0p50(m) && J_0p50(l-1) == JJ_0p50(m)
            TT_0p50(m) = T_0p50(l)-T_0p50(l-1);
            F_0p50(m) = 1/TT_0p50(m);
            JJ_0p50(m) = J_0p50(l);
        else
            m = m+1;
            JJ_0p50(m) = J_0p50(l);
        end
    end
    
    I_0p50 = A*JJ_0p50;

    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_0.75/period.dat
    
    J_0p75 = period(:,1);
    T_0p75 = period(:,2);
    
    N = size(J_0p75,1);
    
    JJ_0p75(1) = J_0p75(1);
    m = 1;
    
    for l = 2:N
        if J_0p75(l) == JJ_0p75(m) && J_0p75(l-1) == JJ_0p75(m)
            TT_0p75(m) = T_0p75(l)-T_0p75(l-1);
            F_0p75(m) = 1/TT_0p75(m);
            JJ_0p75(m) = J_0p75(l);
        else
            m = m+1;
            JJ_0p75(m) = J_0p75(l);
        end
    end
    
    I_0p75 = A*JJ_0p75;

    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_1.00/period.dat
    
    J_1p00 = period(:,1);
    T_1p00 = period(:,2);
    
    N = size(J_1p00,1);
    
    JJ_1p00(1) = J_1p00(1);
    m = 1;
    
    for l = 2:N
        if J_1p00(l) == JJ_1p00(m) && J_1p00(l-1) == JJ_1p00(m)
            TT_1p00(m) = T_1p00(l)-T_1p00(l-1);
            F_1p00(m) = 1/TT_1p00(m);
            JJ_1p00(m) = J_1p00(l);
        else
            m = m+1;
            JJ_1p00(m) = J_1p00(l);
        end
    end
    
    I_1p00 = A*JJ_1p00;

    %----------------------------------------------------------------------
    % \section{2.2 Faz o grafico comparativo entre as frequencias}
    %----------------------------------------------------------------------
    
    figure;
    l_0p25 = plot(I_0p25,F_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,F_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,F_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,F_1p00,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([min([min(F_0p25) min(F_0p50) min(F_0p75) min(F_1p00)]) max([max(F_0p25) max(F_0p50) max(F_0p75) max(F_1p00)])]);
    pos = 'best';
    leg = legend([l_0p25,l_0p50,l_0p75,l_1p00],'$h_s=0.25$','$h_s=0.50$','$h_s=0.75$','$h_s=1.00$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$\nu t_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.035,0.975*max(F_1p00),'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.035,0.94*max(F_1p00),'$l_z=4.8\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('frq','_nt.png');
    file_t = strcat('frq','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0/Lz_4.8'];
    system(cmd);
    

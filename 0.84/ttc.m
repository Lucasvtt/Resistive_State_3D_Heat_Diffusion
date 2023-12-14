%--------------------------------------------------------------------------
% \section{1. Comparacao entre os substratos 
%             para \kappa = 1, T_0 = 0.84T_c e l_z = 3\xi(0)}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    h_s = [0.25 0.50 0.75 1.00]';

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.25/Qm.dat
    
    dy = 0.5;
    dz = 0.5;
    l_y = 20;
    l_z = 3;
    l_y = l_y-dy;
    l_z = l_z-dz;
    A = l_y*l_z;
    
    I_0p25 = A*Qm(:,1);
    Q_0p25 = h_s(1)*Qm(:,2);
    
    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.50/Qm.dat
    
    I_0p50 = A*Qm(:,1);
    Q_0p50 = h_s(2)*Qm(:,2);

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_0.75/Qm.dat
    
    I_0p75 = A*Qm(:,1);
    Q_0p75 = h_s(3)*Qm(:,2);    

    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_1.00/Qm.dat
    
    I_1p00 = A*Qm(:,1);
    Q_1p00 = h_s(4)*Qm(:,2);

    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as 
    %              taxas de tansferencia de calor}
    %----------------------------------------------------------------------
    
    figure;
    l_0p25 = plot(I_0p25,Q_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,Q_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,Q_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,Q_1p00,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([min([min(Q_0p25) min(Q_0p50) min(Q_0p75) min(Q_1p00)]) max([max(Q_0p25) max(Q_0p50) max(Q_0p75) max(Q_1p00)])]);
    pos = 'southwest';
    leg = legend([l_0p25,l_0p50,l_0p75,l_1p00],'$h_s=0.25$','$h_s=0.50$','$h_s=0.75$','$h_s=1.00$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$\dot{Q}/\dot{Q}_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.1,0.9*max(Q_1p00),'$T_0=0.84T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.1,0.8*max(Q_1p00),'$l_z=3\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('Q','_nt.png');
    file_t = strcat('Q','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.84/kappa_1.0/Lz_3'];
    system(cmd);

%--------------------------------------------------------------------------
% \section{2. Comparacao entre os substratos 
%             para \kappa = 1, T_0 = 0.84T_c e l_z = 6\xi(0)}
%--------------------------------------------------------------------------

    clear all;
%     close all;
    
    %----------------------------------------------------------------------
    % \subsection{2.1 Carrega os arquivos}
    %----------------------------------------------------------------------

    h_s = [0.25 0.50 0.75 1.00]';
    
    load T_0.84/kappa_1.0/Lz_6/hf_0.25_hs_0.25/Qm.dat
    
    dy = 0.5;
    dz = 0.5;
    l_y = 20;
    l_z = 6;
    l_y = l_y-dy;
    l_z = l_z-dz;
    A = l_y*l_z;
    
    I_0p25 = A*Qm(:,1);
    Q_0p25 = h_s(1)*Qm(:,2);
    
    load T_0.84/kappa_1.0/Lz_6/hf_0.25_hs_0.50/Qm.dat
    
    I_0p50 = A*Qm(:,1);
    Q_0p50 = h_s(2)*Qm(:,2);

    load T_0.84/kappa_1.0/Lz_6/hf_0.25_hs_0.75/Qm.dat
    
    I_0p75 = A*Qm(:,1);
    Q_0p75 = h_s(3)*Qm(:,2);    

    load T_0.84/kappa_1.0/Lz_6/hf_0.25_hs_1.00/Qm.dat
    
    I_1p00 = A*Qm(:,1);
    Q_1p00 = h_s(4)*Qm(:,2);

    %----------------------------------------------------------------------
    % \section{2.2 Faz o grafico comparativo entre as 
    %              taxas de tansferencia de calor}
    %----------------------------------------------------------------------
    
    figure;
    l_0p25 = plot(I_0p25,Q_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,Q_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,Q_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,Q_1p00,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([min([min(Q_0p25) min(Q_0p50) min(Q_0p75) min(Q_1p00)]) max([max(Q_0p25) max(Q_0p50) max(Q_0p75) max(Q_1p00)])]);
    pos = 'southwest';
    leg = legend([l_0p25,l_0p50,l_0p75,l_1p00],'$h_s=0.25$','$h_s=0.50$','$h_s=0.75$','$h_s=1.00$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$\dot{Q}/\dot{Q}_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.15,0.9*max(Q_1p00),'$T_0=0.84T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.15,0.8*max(Q_1p00),'$l_z=6\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('Q','_nt.png');
    file_t = strcat('Q','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.84/kappa_1.0/Lz_6'];
    system(cmd);    
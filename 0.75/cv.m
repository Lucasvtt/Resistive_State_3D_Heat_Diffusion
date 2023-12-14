%--------------------------------------------------------------------------
% \section{1. Comparacao entre os substratos 
%             para \kappa = 1, T_0 = 0.75T_c e l_z = 2.4\xi(0)}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.25/IV.dat
    
    I_0p25 = IV(:,2);
    V_0p25 = IV(:,3);
    R_0p25 = diffxy(I_0p25,V_0p25);
    fct = 3/max(R_0p25);
    
    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.50/IV.dat
    
    I_0p50 = IV(:,2);
    V_0p50 = IV(:,3);
    R_0p50 = diffxy(I_0p50,V_0p50);

    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.75/IV.dat
    
    I_0p75 = IV(:,2);
    V_0p75 = IV(:,3);
    R_0p75 = diffxy(I_0p75,V_0p75);    

    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/IV.dat
    
    I_1p00 = IV(:,2);
    V_1p00 = IV(:,3);    
    R_1p00 = diffxy(I_1p00,V_1p00);

    dy = 0.4;
    dz = 0.4;
    l_x = 24;
    l_y = 16;
    l_z = 2.4;
    A = (l_y-dy)*(l_z-dz);
    V_o = l_x*IV(:,1);  % voltagem Ohmica
    sigma = 1;
    R_o = l_x/(A*sigma)*I_1p00./I_1p00;
    V_o = R_o.*I_1p00;
    
    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_0p25 = plot(I_0p25,V_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,V_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,V_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,V_1p00,'LineWidth',1.5);
    l_o = plot(I_1p00,V_o,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([min([min(V_0p25) min(V_0p50) min(V_0p75) min(V_1p00)]) max([max(V_0p25) max(V_0p50) max(V_0p75) max(V_1p00)])]);
    pos = 'west';
    leg = legend([l_0p25,l_0p50,l_0p75,l_1p00,l_o],'$h_s=0.25$','$h_s=0.50$','$h_s=0.75$','$h_s=1.00$','Ohm','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    set(gca, 'XTickLabel', [' ' '0.5' '1.0' '1.5' '2.0']);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.1,0.98*max(V_0p25),'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.1,0.88*max(V_0p25),'$l_z=2.4\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV','_nt.png');
    file_t = strcat('IV','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0/Lz_2.4'];
    system(cmd);

    %----------------------------------------------------------------------
    % \section{1.3 Faz o grafico comparativo entre as resistencias}
    %----------------------------------------------------------------------

    figure;
    l_0p25 = plot(I_0p25,R_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,R_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,R_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,R_1p00,'LineWidth',1.5);
%     l_o = plot(I_1p00,R_o,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([0 fct*max([max(R_0p25) fct*max(R_0p50) fct*max(R_0p75) fct*max(R_1p00)])]);
    pos = 'west';
    leg = legend([l_0p25,l_0p50,l_0p75,l_1p00],'$h_s=0.25$','$h_s=0.50$','$h_s=0.75$','$h_s=1.00$','Location',pos);    
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$(dV/dI)/(\varphi_{GL}/I_{GL})$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.1,0.85*fct*max(R_0p25),'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.1,0.75*fct*max(R_0p25),'$l_z=2.4\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IR','_nt.png');
    file_t = strcat('IR','_t.png');
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
    
    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_0.25/IV.dat
    
    I_0p25 = IV(:,2);
    V_0p25 = IV(:,3);
    R_0p25 = diffxy(I_0p25,V_0p25);
    fct = 2/max(R_0p25);
    
    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_0.50/IV.dat
    
    I_0p50 = IV(:,2);
    V_0p50 = IV(:,3);
    R_0p50 = diffxy(I_0p50,V_0p50);

    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_0.75/IV.dat
    
    I_0p75 = IV(:,2);
    V_0p75 = IV(:,3);
    R_0p75 = diffxy(I_0p75,V_0p75);    

    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_1.00/IV.dat
    
    I_1p00 = IV(:,2);
    V_1p00 = IV(:,3);    
    R_1p00 = diffxy(I_1p00,V_1p00);
    fct = 1.5/max(R_1p00);

    dy = 0.4;
    dz = 0.4;
    l_x = 24;
    l_y = 16;
    l_z = 4.8;
    A = (l_y-dy)*(l_z-dz);
    V_o = l_x*IV(:,1);  % voltagem Ohmica
    sigma = 1;
    R_o = l_x/(A*sigma)*I_1p00./I_1p00;
    V_o = R_o.*I_1p00;
    
    %----------------------------------------------------------------------
    % \section{2.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_0p25 = plot(I_0p25,V_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,V_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,V_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,V_1p00,'LineWidth',1.5);
    l_o = plot(I_1p00,V_o,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([min([min(V_0p25) min(V_0p50) min(V_0p75) min(V_1p00)]) max([max(V_0p25) max(V_0p50) max(V_0p75) max(V_1p00)])]);
    pos = 'west';
    leg = legend([l_0p25,l_0p50,l_0p75,l_1p00,l_o],'$h_s=0.25$','$h_s=0.50$','$h_s=0.75$','$h_s=1.00$','Ohm','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.2,0.9*max(V_0p75),'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.2,0.8*max(V_0p75),'$l_z=4.8\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV','_nt.png');
    file_t = strcat('IV','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0/Lz_4.8'];
    system(cmd);

    %----------------------------------------------------------------------
    % \section{2.3 Faz o grafico comparativo entre as resistencias}
    %----------------------------------------------------------------------

    figure;
    l_0p25 = plot(I_0p25,R_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,R_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,R_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,R_1p00,'LineWidth',1.5);
%     l_o = plot(I_1p00,R_o,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([0 fct*max([max(R_0p25) fct*max(R_0p50) fct*max(R_0p75) fct*max(R_1p00)])]);
    pos = 'west';
    leg = legend([l_0p25,l_0p50,l_0p75,l_1p00],'$h_s=0.25$','$h_s=0.50$','$h_s=0.75$','$h_s=1.00$','Location',pos);    
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$(dV/dI)/(\varphi_{GL}/I_{GL})$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.2,0.95*max(R_1p00)/6,'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.2,0.85*max(R_1p00)/6,'$l_z=4.8\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IR','_nt.png');
    file_t = strcat('IR','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0/Lz_4.8'];
    system(cmd);   
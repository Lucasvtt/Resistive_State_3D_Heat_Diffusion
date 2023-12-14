%--------------------------------------------------------------------------
% \section{1. Comparacao entre as espessuras 
%             para \kappa = 1, T_0 = 0.75T_c, h_f = 0.25, h_s = 1.00}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/IV.dat
    
    I_2p4 = IV(:,2);
    V_2p4 = IV(:,3);    
    R_2p4 = diffxy(I_2p4,V_2p4);
    fct = 0.5;

    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_1.00/IV.dat
    
    I_4p8 = IV(:,2);
    V_4p8 = IV(:,3);    
    R_4p8 = diffxy(I_4p8,V_4p8);    
    
    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_2p4 = plot(I_2p4,V_2p4,'LineWidth',1.5);
    hold on;
    l_4p8 = plot(I_4p8,V_4p8,'LineWidth',1.5);
    xlim([min([min(I_2p4) min(I_4p8)]) max([max(I_2p4) max(I_4p8)])]);
    ylim([min([min(V_2p4) min(V_4p8)]) max([max(V_2p4) max(V_4p8)])]);
    pos = 'northwest';
    leg = legend([l_2p4,l_4p8],'$l_z=2.4\xi(0)$','$l_z=4.8\xi(0)$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_2p4)+0.25,0.5*max(V_2p4),'$T_0=0.75T_c$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_2p4)+0.25,0.4*max(V_2p4),'$\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_2p4)+0.25,0.3*max(V_2p4),'$h_s=1$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV_thick','_nt.png');
    file_t = strcat('IV_thick','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0'];
    system(cmd);

    %----------------------------------------------------------------------
    % \section{1.3 Faz o grafico comparativo entre as resistencias}
    %----------------------------------------------------------------------

    figure;
    l_2p4 = plot(I_2p4,R_2p4,'LineWidth',1.5);
    hold on;
    l_4p8 = plot(I_4p8,R_4p8,'LineWidth',1.5);
    xlim([min([min(I_2p4) min(I_4p8)]) max([max(I_2p4) max(I_4p8)])]);
    ylim(fct*[min([min(R_2p4) min(R_4p8)]) fct*max([max(R_2p4) max(R_4p8)])]);
    pos = 'northwest';
    leg = legend([l_2p4,l_4p8],'$l_z=2.4\xi(0)$','$l_z=4.8\xi(0)$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$(dV/dI)/(\varphi_{GL}/I_{GL})$','Interpreter','LaTeX','FontSize',14);
    text(min(I_2p4)+0.25,0.25*fct*max(R_2p4),'$T_0=0.75T_c$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_2p4)+0.25,0.2*fct*max(R_2p4),'$\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_2p4)+0.25,0.15*fct*max(R_2p4),'$h_s=1$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IR_thick','_nt.png');
    file_t = strcat('IR_thick','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0'];
    system(cmd);        
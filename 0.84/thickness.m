%--------------------------------------------------------------------------
% \section{1. Comparacao entre as espessuras 
%             para \kappa = 1, T_0 = 0.84T_c, h_f = 0.25, h_s = 1.00}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_1.00/IV.dat
    
    I_3 = IV(:,2);
    V_3 = IV(:,3);    
    R_3 = diffxy(I_3,V_3);
    fct = 0.5;

    load T_0.84/kappa_1.0/Lz_6/hf_0.25_hs_1.00/IV.dat
    
    I_6 = IV(:,2);
    V_6 = IV(:,3);    
    R_6 = diffxy(I_6,V_6);    
    
    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_3 = plot(I_3,V_3,'LineWidth',1.5);
    hold on;
    l_6 = plot(I_6,V_6,'LineWidth',1.5);
    xlim([min([min(I_3) min(I_6)]) max([max(I_3) max(I_6)])]);
    ylim([min([min(V_3) min(V_6)]) max([max(V_3) max(V_6)])]);
    pos = 'northwest';
    leg = legend([l_3,l_6],'$l_z=3\xi(0)$','$l_z=6\xi(0)$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_3)+0.25,0.5*max(V_3),'$T_0=0.84T_c$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_3)+0.25,0.4*max(V_3),'$\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_3)+0.25,0.3*max(V_3),'$h_s=1$','FontSize',24,'Interpreter','LaTeX');
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
    
    cmd = ['mv -f *.png', ' ', 'T_0.84/kappa_1.0'];
    system(cmd);

    %----------------------------------------------------------------------
    % \section{1.3 Faz o grafico comparativo entre as resistencias}
    %----------------------------------------------------------------------

    figure;
    l_3 = plot(I_3,R_3,'LineWidth',1.5);
    hold on;
    l_6 = plot(I_6,R_6,'LineWidth',1.5);
    xlim([min([min(I_3) min(I_6)]) max([max(I_3) max(I_6)])]);
    ylim(fct*[min([min(R_3) min(R_6)]) fct*max([max(R_3) max(R_6)])]);
    pos = 'northwest';
    leg = legend([l_3,l_6],'$l_z=3\xi(0)$','$l_z=6\xi(0)$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$(dV/dI)/(\varphi_{GL}/I_{GL})$','Interpreter','LaTeX','FontSize',14);
    text(min(I_3)+0.25,0.25*fct*max(R_3),'$T_0=0.84T_c$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_3)+0.25,0.2*fct*max(R_3),'$\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_3)+0.25,0.15*fct*max(R_3),'$h_s=1$','FontSize',24,'Interpreter','LaTeX');
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
    
    cmd = ['mv -f *.png', ' ', 'T_0.84/kappa_1.0'];
    system(cmd);
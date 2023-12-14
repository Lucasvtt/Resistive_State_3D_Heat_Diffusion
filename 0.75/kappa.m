%--------------------------------------------------------------------------
% \section{1. Comparacao entre os kappa's 
%             para \kappa = 1, \kappa = 2, T_0 = 0.75T_c e l_z = 2.4\xi(0)}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
  
    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/IV.dat
    
    I_k_1 = IV(:,2);
    V_k_1 = IV(:,3);    
    R_k_1 = diffxy(I_k_1,V_k_1);

    load T_0.75/kappa_2.0/Lz_2.4/hf_0.25_hs_1.00/IV.dat
    
    I_k_2 = IV(:,2);
    V_k_2 = IV(:,3);    
    R_k_2 = diffxy(I_k_2,V_k_2);
    fct = 0.3;
    
    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_k_1 = plot(I_k_1,V_k_1,'LineWidth',1.5);
    hold on;
    l_k_2 = plot(I_k_2,V_k_2,'LineWidth',1.5);
    xlim([min([min(I_k_1) min(I_k_2)]) max([max(I_k_1) max(I_k_2)])]);
    ylim([min([min(V_k_1) min(V_k_2)]) max([max(V_k_1) max(V_k_2)])]);
    pos = 'northwest';
    leg = legend([l_k_1,l_k_2],'$\kappa=1$','$\kappa=2$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_k_1)+0.1,0.6*max(V_k_1),'$T_0=0.75T_c,\,h_s=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_k_1)+0.1,0.5*max(V_k_1),'$l_z=2.4\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV_kappa_l_z_2p4','_nt.png');
    file_t = strcat('IV_kappa_l_z_2p4','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75'];
    system(cmd);

    %----------------------------------------------------------------------
    % \section{1.3 Faz o grafico comparativo entre as resistencias}
    %----------------------------------------------------------------------

    figure;
    l_k_1 = plot(I_k_1,R_k_1,'LineWidth',1.5);
    hold on;
    l_k_2 = plot(I_k_2,R_k_2,'LineWidth',1.5);
    xlim([min([min(I_k_1) min(I_k_2)]) max([max(I_k_1) max(I_k_2)])]);
    ylim(fct*[min([min(R_k_1) min(R_k_2)]) max([max(R_k_1) max(R_k_2)])]);
    pos = 'northwest';
    leg = legend([l_k_1,l_k_2],'$\kappa=1$','$\kappa=2$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$(dV/dI)(\varphi_{GL}/I_{GL})$','Interpreter','LaTeX','FontSize',14);
    text(min(I_k_1)+0.1,0.6*fct*max(R_k_2),'$T_0=0.75T_c,\,h_s=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_k_1)+0.1,0.5*fct*max(R_k_2),'$l_z=2.4\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IR_kappa_l_z_2p4','_nt.png');
    file_t = strcat('IR_kappa_l_z_2p4','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75'];
    system(cmd);

%--------------------------------------------------------------------------
% \section{2. Comparacao entre os kappa's 
%             para \kappa = 1, \kappa = 2, T_0 = 0.75T_c e l_z = 4.8\xi(0)}
%--------------------------------------------------------------------------

    clear all;
%     close all;
    
    %----------------------------------------------------------------------
    % \subsection{2.1 Carrega os arquivos}
    %----------------------------------------------------------------------
  
    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_1.00/IV.dat
    
    I_k_1 = IV(:,2);
    V_k_1 = IV(:,3);    
    R_k_1 = diffxy(I_k_1,V_k_1);

    load T_0.75/kappa_2.0/Lz_4.8/hf_0.25_hs_1.00/IV.dat
    
    I_k_2 = IV(:,2);
    V_k_2 = IV(:,3);    
    R_k_2 = diffxy(I_k_2,V_k_2);
    fct = 0.3;
    
    %----------------------------------------------------------------------
    % \section{2.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_k_1 = plot(I_k_1,V_k_1,'LineWidth',1.5);
    hold on;
    l_k_2 = plot(I_k_2,V_k_2,'LineWidth',1.5);
    xlim([min([min(I_k_1) min(I_k_2)]) max([max(I_k_1) max(I_k_2)])]);
    ylim([min([min(V_k_1) min(V_k_2)]) max([max(V_k_1) max(V_k_2)])]);
    pos = 'northwest';
    leg = legend([l_k_1,l_k_2],'$\kappa=1$','$\kappa=2$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_k_1)+0.1,0.6*max(V_k_1),'$T_0=0.75T_c,\,h_s=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_k_1)+0.1,0.5*max(V_k_1),'$l_z=4.8\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV_kappa_l_z_4p8','_nt.png');
    file_t = strcat('IV_kappa_l_z_4p8','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75'];
    system(cmd);

    %----------------------------------------------------------------------
    % \section{1.3 Faz o grafico comparativo entre as resistencias}
    %----------------------------------------------------------------------

    figure;
    l_k_1 = plot(I_k_1,R_k_1,'LineWidth',1.5);
    hold on;
    l_k_2 = plot(I_k_2,R_k_2,'LineWidth',1.5);
    xlim([min([min(I_k_1) min(I_k_2)]) max([max(I_k_1) max(I_k_2)])]);
    ylim(fct*[min([min(R_k_1) min(R_k_2)]) max([max(R_k_1) max(R_k_2)])]);
    pos = 'northwest';
    leg = legend([l_k_1,l_k_2],'$\kappa=1$','$\kappa=2$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$(dV/dI)(\varphi_{GL}/I_{GL})$','Interpreter','LaTeX','FontSize',14);
    text(min(I_k_1)+0.1,0.6*fct*max(R_k_2),'$T_0=0.75T_c,\,h_s=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_k_1)+0.1,0.5*fct*max(R_k_2),'$l_z=4.8\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IR_kappa_l_z_4p8','_nt.png');
    file_t = strcat('IR_kappa_l_z_4p8','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75'];
    system(cmd);     

%--------------------------------------------------------------------------
% \section{3. Comparacao entre os kappa's 
%             para \kappa = 1, \kappa = 2, T_0 = 0.75T_c e l_z = 2.4\xi(0)}
%--------------------------------------------------------------------------    

    clear all;
%     close all;
    
    %----------------------------------------------------------------------
    % \subsection{3.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/period.dat
    
    dy = 0.4;
    dz = 0.4;
    l_y = 16;
    l_z = 2.4;
    l_y = l_y-dy;
    l_z = l_z-dz;
    A = l_y*l_z;
    
    J_k_1 = period(:,1); 
    T_k_1 = period(:,2);
      
    N = size(J_k_1,1);
    
    JJ_k_1(1) = J_k_1(1);
    m = 1;
    
    for l = 2:N
        if J_k_1(l) == JJ_k_1(m) && J_k_1(l-1) == JJ_k_1(m)
            TT_k_1(m) = T_k_1(l)-T_k_1(l-1);
            F_k_1(m) = 1/TT_k_1(m);
            JJ_k_1(m) = J_k_1(l);
        else
            m = m+1;
            JJ_k_1(m) = J_k_1(l);
        end
    end
    
    I_k_1 = A*JJ_k_1;
    
    load T_0.75/kappa_2.0/Lz_2.4/hf_0.25_hs_1.00/period.dat
    
    J_k_2 = period(:,1); 
    T_k_2 = period(:,2);
      
    N = size(J_k_2,1);
    
    JJ_k_2(1) = J_k_2(1);
    m = 1;
    
    for l = 2:N
        if J_k_2(l) == JJ_k_2(m) && J_k_2(l-1) == JJ_k_2(m)
            TT_k_2(m) = T_k_2(l)-T_k_2(l-1);
            F_k_2(m) = 1/TT_k_2(m);
            JJ_k_2(m) = J_k_2(l);
        else
            m = m+1;
            JJ_k_2(m) = J_k_2(l);
        end
    end
    
    I_k_2 = A*JJ_k_2;

    %----------------------------------------------------------------------
    % \section{3.2 Faz o grafico comparativo entre as frequencias}
    %----------------------------------------------------------------------
    
    figure;
    l_k_1 = plot(I_k_1,F_k_1,'LineWidth',1.5);
    hold on;
    l_k_2 = plot(I_k_2,F_k_2,'LineWidth',1.5);
    xlim([min([min(I_k_1) min(I_k_2)]) max([max(I_k_1) max(I_k_2)])]);
    ylim([min([min(F_k_1) min(F_k_2)]) max([max(F_k_1) max(F_k_2)])]);
    pos = 'northwest';
    leg = legend([l_k_1,l_k_2],'$\kappa=1$','$\kappa=2$','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$\nu t_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_k_1)+0.3,min(F_k_1)+0.035,'$T_0=0.75T_c,\,h_s=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_k_1)+0.3,min(F_k_1)+0.02,'$l_z=2.4\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('frq_kappa','_nt.png');
    file_t = strcat('frq_kappa','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75'];
    system(cmd);    
%--------------------------------------------------------------------------
% \section{1. Comparacao entre substrato forte e sem difusao de calor 
%             para \kappa = 1, T_0 = 0.84T_c e l_z = 3\xi(0)}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.84/kappa_1.0/Lz_3/hf_0.25_hs_1.00/IV.dat
    
    I_1p00 = IV(:,2);
    V_1p00 = IV(:,3);    
    R_1p00 = diffxy(I_1p00,V_1p00);
    fct = 2/max(R_1p00);

    load T_0.84/kappa_1.0/Lz_3/ideal/IV.dat
    
    I_id = IV(:,2);
    V_id = IV(:,3);    
    R_id = diffxy(I_id,V_id);  

    load T_0.84/kappa_1.0/Lz_3/No_Heat_Diffusion/IV.dat
    
    I_nhd = IV(:,2);
    V_nhd = IV(:,3);    
    R_nhd = diffxy(I_nhd,V_nhd); 

    dy = 0.5;
    dz = 0.5;
    l_x = 30;
    l_y = 20;
    l_z = 3;
    A = (l_y-dy)*(l_z-dz);
    V_o = l_x*IV(:,1);  % voltagem Ohmica
    sigma = 1;
    R_o = l_x/(A*sigma)*I_1p00./I_1p00;  
    
    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_1p00 = plot(I_1p00,V_1p00,'LineWidth',1.5);
    hold on; 
    l_id = plot(I_id,V_id,'LineWidth',1.5); 
    l_nhd = plot(I_nhd,V_nhd,'LineWidth',1.5);
    l_o = plot(I_nhd,V_o,'LineWidth',1.5);
    xlim([min([min(I_1p00) min(I_id) min(I_nhd)]) max([max(I_1p00) max(I_id) max(I_nhd)])]);
    ylim([min([min(V_1p00) min(V_id) min(V_nhd)]) max([max(V_1p00) max(V_id) max(V_nhd)])]);
    pos = 'northwest';
    leg = legend([l_1p00,l_id,l_nhd,l_o],'$h_s=1.00$','Ideal','No heat diffusion','Ohm','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_1p00)+0.1,0.55*max(V_nhd),'$T_0=0.84T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_1p00)+0.1,0.45*max(V_nhd),'$l_z=3\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV_nhd','_nt.png');
    file_t = strcat('IV_nhd','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.84/kappa_1.0/Lz_3'];
    system(cmd);

    %----------------------------------------------------------------------
    % \section{1.3 Faz o grafico comparativo entre as resistencias}
    %----------------------------------------------------------------------

    figure;
    l_1p00 = plot(I_1p00,R_1p00,'LineWidth',1.5);
    hold on;
    l_id = plot(I_id,R_id,'LineWidth',1.5); 
    l_nhd = plot(I_nhd,R_nhd,'LineWidth',1.5);
    xlim([min([min(I_1p00) min(I_id) min(I_nhd)]) max([max(I_1p00) max(I_id) max(I_nhd)])]);
    ylim([fct*min([min(R_1p00) min(R_id) min(R_nhd)]) fct*max([max(R_1p00) max(R_id) max(R_nhd)])]);
    pos = 'northwest';
    leg = legend([l_1p00,l_id,l_nhd],'$h_s=1.00$','Ideal','No heat diffusion','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$(dV/dI)/(\varphi_{GL}/I_{GL})$','Interpreter','LaTeX','FontSize',14);
    text(min(I_1p00)+0.1,0.55*fct*max(R_nhd),'$T_0=0.84T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_1p00)+0.1,0.45*fct*max(R_nhd),'$l_z=3\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_1p00)+0.2,0.25*fct*max(R_nhd),'Meissner State','FontSize',24,'Interpreter','LaTeX');
    text(min(I_1p00)+1.3,0.25*fct*max(R_nhd),'Resistive State','FontSize',24,'Interpreter','LaTeX');    
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IR_nhd','_nt.png');
    file_t = strcat('IR_nhd','_t.png');
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
% \section{2. Comparacao entre substrato forte e sem difusao de calor 
%             para \kappa = 1, T_0 = 0.84T_c e l_z = 6\xi(0)}
%--------------------------------------------------------------------------

    clear all;
%     close all;
    
    %----------------------------------------------------------------------
    % \subsection{2.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.84/kappa_1.0/Lz_6/hf_0.25_hs_1.00/IV.dat
    
    I_1p00 = IV(:,2);
    V_1p00 = IV(:,3);    
    R_1p00 = diffxy(I_1p00,V_1p00);
    fct = 2/max(R_1p00);

    load T_0.84/kappa_1.0/Lz_6/ideal/IV.dat
    
    I_id = IV(:,2);
    V_id = IV(:,3);    
    R_id = diffxy(I_id,V_id);  

    load T_0.84/kappa_1.0/Lz_6/No_Heat_Diffusion/IV.dat
    
    I_nhd = IV(:,2);
    V_nhd = IV(:,3);    
    R_nhd = diffxy(I_nhd,V_nhd); 

    dy = 0.5;
    dz = 0.5;
    l_x = 30;
    l_y = 20;
    l_z = 6;
    A = (l_y-dy)*(l_z-dz);
    V_o = l_x*IV(:,1);  % voltagem Ohmica
    sigma = 1;
    R_o = l_x/(A*sigma)*I_1p00./I_1p00;  
    
    %----------------------------------------------------------------------
    % \section{2.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_1p00 = plot(I_1p00,V_1p00,'LineWidth',1.5);
    hold on; 
    l_id = plot(I_id,V_id,'LineWidth',1.5); 
    l_nhd = plot(I_nhd,V_nhd,'LineWidth',1.5);
    l_o = plot(I_nhd,V_o,'LineWidth',1.5);
    xlim([min([min(I_1p00) min(I_id) min(I_nhd)]) max([max(I_1p00) max(I_id) max(I_nhd)])]);
    ylim([min([min(V_1p00) min(V_id) min(V_nhd)]) max([max(V_1p00) max(V_nhd) max(V_id)])]);
    pos = 'northwest';
    leg = legend([l_1p00,l_id,l_nhd,l_o],'$h_s=1.00$','Ideal','No heat diffusion','Ohm','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_1p00)+0.15,0.55*max(V_nhd),'$T_0=0.84T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_1p00)+0.15,0.45*max(V_nhd),'$l_z=6\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV_nhd','_nt.png');
    file_t = strcat('IV_nhd','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.84/kappa_1.0/Lz_6'];
    system(cmd);

    %----------------------------------------------------------------------
    % \section{2.3 Faz o grafico comparativo entre as resistencias}
    %----------------------------------------------------------------------

    figure;
    l_1p00 = plot(I_1p00,R_1p00,'LineWidth',1.5);
    hold on;
    l_id = plot(I_id,R_id,'LineWidth',1.5); 
    l_nhd = plot(I_nhd,R_nhd,'LineWidth',1.5);
    xlim([min([min(I_1p00) min(I_id) min(I_nhd)]) max([max(I_1p00) max(I_id) max(I_nhd)])]);
    ylim([fct*min([min(R_1p00) min(R_id) min(R_nhd)]) fct*max([max(R_1p00) max(R_id) max(R_nhd)])]);
    pos = 'northwest';
    leg = legend([l_1p00,l_id,l_nhd],'$h_s=1.00$','Ideal','No heat diffusion','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$(dV/dI)/(\varphi_{GL}/I_{GL})$','Interpreter','LaTeX','FontSize',14);
    text(min(I_1p00)+0.2,0.71*fct*max(R_nhd),'$T_0=0.84T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_1p00)+0.2,0.58*fct*max(R_nhd),'$l_z=6\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IR_nhd','_nt.png');
    file_t = strcat('IR_nhd','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.84/kappa_1.0/Lz_6'];
    system(cmd);    
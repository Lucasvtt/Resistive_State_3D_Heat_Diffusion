%--------------------------------------------------------------------------
% \section{1. Comparacao entre os substratos (curvas de histere) 
%             para \kappa = 1, T_0 = 0.75T_c e l_z = 2.4\xi(0)}
%--------------------------------------------------------------------------

    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/IV.dat
    
    I_1p00 = IV(:,2);
    V_1p00 = IV(:,3);    

    load T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/Down_Sweep/IV.dat
    
    I_1p00_ds = IV(:,2);
    V_1p00_ds = IV(:,3);

    load T_0.75/kappa_1.0/Lz_2.4/ideal/IV.dat
    
    I_id = IV(:,2);
    V_id = IV(:,3);    
    
    load T_0.75/kappa_1.0/Lz_2.4/ideal/Down_Sweep/IV.dat
    
    I_id_ds = IV(:,2);
    V_id_ds = IV(:,3);    

    load T_0.75/kappa_1.0/Lz_2.4/No_Heat_Diffusion/IV.dat
    
    I_nhd = IV(:,2);
    V_nhd = IV(:,3); 

    load T_0.75/kappa_1.0/Lz_2.4/No_Heat_Diffusion/Down_Sweep/IV.dat
    
    I_nhd_ds = IV(:,2);
    V_nhd_ds = IV(:,3);    
    
    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_1p00 = plot(I_1p00,V_1p00,'Color',[0 0.4470 0.7410],'LineWidth',1.5);
    hold on;
    plot(I_1p00_ds,V_1p00_ds,'--','Color',[0 0.4470 0.7410],'LineWidth',1.5);
    l_id = plot(I_id,V_id,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
    plot(I_id_ds,V_id_ds,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
    l_nhd = plot(I_nhd,V_nhd,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
    plot(I_nhd_ds,V_nhd_ds,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);    
    xlim([min([min(I_1p00) min(I_id) min(I_nhd)]) max([max(I_1p00) max(I_id) max(I_nhd)])]);
    ylim([min([min(V_1p00) min(V_id) min(V_nhd)]) max([max(V_1p00) max(V_id) max(V_nhd)])]);
    pos = 'northwest';
    leg = legend([l_1p00 l_id,l_nhd],'$h_s=1.00$','ideal','No heat diffusion','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_1p00)+0.2,0.6*max(V_nhd),'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_1p00)+0.2,0.5*max(V_nhd),'$l_z=2.4\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV_histeresis','_nt.png');
    file_t = strcat('IV_histeresis','_t.png');
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
% \section{2. Comparacao entre os substratos (curvas de histere) 
%             para \kappa = 1, T_0 = 0.75T_c e l_z = 6\xi(0)}
%--------------------------------------------------------------------------

    clear all;
%     close all;
    
    %----------------------------------------------------------------------
    % \subsection{2.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_1.00/IV.dat
    
    I_1p00 = IV(:,2);
    V_1p00 = IV(:,3);    

    load T_0.75/kappa_1.0/Lz_4.8/hf_0.25_hs_1.00/Down_Sweep/IV.dat
    
    I_1p00_ds = IV(:,2);
    V_1p00_ds = IV(:,3);

    load T_0.75/kappa_1.0/Lz_4.8/ideal/IV.dat
    
    I_id = IV(:,2);
    V_id = IV(:,3);    
    
    load T_0.75/kappa_1.0/Lz_4.8/ideal/Down_Sweep/IV.dat
    
    I_id_ds = IV(:,2);
    V_id_ds = IV(:,3);    

    load T_0.75/kappa_1.0/Lz_4.8/No_Heat_Diffusion/IV.dat
    
    I_nhd = IV(:,2);
    V_nhd = IV(:,3); 

    load T_0.75/kappa_1.0/Lz_4.8/No_Heat_Diffusion/Down_Sweep/IV.dat
    
    I_nhd_ds = IV(:,2);
    V_nhd_ds = IV(:,3);    
    
    %----------------------------------------------------------------------
    % \section{2.2 Faz o grafico comparativo entre as voltagens}
    %----------------------------------------------------------------------
    
    figure;
    l_1p00 = plot(I_1p00,V_1p00,'Color',[0 0.4470 0.7410],'LineWidth',1.5);
    hold on;
    plot(I_1p00_ds,V_1p00_ds,'--','Color',[0 0.4470 0.7410],'LineWidth',1.5);
    l_id = plot(I_id,V_id,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
    plot(I_id_ds,V_id_ds,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
    l_nhd = plot(I_nhd,V_nhd,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
    plot(I_nhd_ds,V_nhd_ds,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);    
    xlim([min([min(I_1p00) min(I_id) min(I_nhd)]) max([max(I_1p00) max(I_id) max(I_nhd)])]);
    ylim([min([min(V_1p00) min(V_id) min(V_nhd)]) max([max(V_1p00) max(V_id) max(V_nhd)])]);
    pos = 'northwest';
    leg = legend([l_1p00 l_id,l_nhd],'$h_s=1.00$','ideal','No heat diffusion','Location',pos);
    legend('boxoff');
    set(leg,'Interpreter','LaTeX');
    %set(gca, 'XTickLabel', []);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'TickLength',[0.02, 0.01]);
    set(gca,'LineWidth',1.5);
    set(gca,'TickDir','in');
    xlabel('$I/I_{GL}$','Interpreter','LaTeX','FontSize',14);
    ylabel('$V/\varphi_{GL}$','Interpreter','LaTeX','FontSize',14);
    text(min(I_1p00)+0.2,0.6*max(V_nhd),'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_1p00)+0.2,0.5*max(V_nhd),'$l_z=4.8\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('IV_histeresis','_nt.png');
    file_t = strcat('IV_histeresis','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0/Lz_4.8'];
    system(cmd);    

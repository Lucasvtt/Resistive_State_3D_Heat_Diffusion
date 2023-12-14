%--------------------------------------------------------------------------
% \section{1. Comparacao entre as temperaturas
%             para \kappa = 1, T_0 = 0.75T_c e l_z = 2.4\xi(0)}
%--------------------------------------------------------------------------


    clear all;
    close all;
    
    %----------------------------------------------------------------------
    % \subsection{1.1 Carrega os arquivos}
    %----------------------------------------------------------------------
    
    dy = 0.4;
    dz = 0.4;
    l_y = 16;
    l_z = 2.4;
    l_y = l_y-dy;
    l_z = l_z-dz;
    A = l_y*l_z;

    T_0 = 0.75;

    dJ_a = 0.0005;

    J_a = 0;

    J_a_max = 0.0615;

    l = 0;

    Dir = 'T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.25/Tx/';
    
    while J_a < J_a_max

        l = l+1;
        J_a_string = num2str(J_a);
    
        if length(J_a_string) == 7
            eval(['load ',Dir,J_a_string,'/ttx.dat']);
        elseif length(J_a_string) == 6
            eval(['load ',Dir,J_a_string,num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 5
            eval(['load ',Dir,J_a_string,num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 4
            eval(['load ',Dir,J_a_string,num2str(0),num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 3
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 2
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);            
        elseif length(J_a_string) == 1
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);            
        else
            eval(['load ',Dir,J_a_string,'/ttx.dat']);
        end
    
        Temp_max(l) = max(ttx(:,2));
        Temp_max(l) = T_0+Temp_max(l); % temperatura maxima no centro do vortice fechado
        curr(l,1) = J_a; 
    
        J_a = J_a+dJ_a;

    end
    
    I_0p25 = A*curr';
    T_0p25 = Temp_max';

    J_a = 0;

    J_a_max = 0.066;

    l = 0;    

    Dir = 'T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.50/Tx/';
    
    while J_a < J_a_max

        l = l+1;
        J_a_string = num2str(J_a);
    
        if length(J_a_string) == 7
            eval(['load ',Dir,J_a_string,'/ttx.dat']);
        elseif length(J_a_string) == 6
            eval(['load ',Dir,J_a_string,num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 5
            eval(['load ',Dir,J_a_string,num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 4
            eval(['load ',Dir,J_a_string,num2str(0),num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 3
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 2
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);            
        elseif length(J_a_string) == 1
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);            
        else
            eval(['load ',Dir,J_a_string,'/ttx.dat']);
        end
    
        Temp_max(l) = max(ttx(:,2));
        Temp_max(l) = T_0+Temp_max(l); % temperatura maxima no centro do vortice fechado
        curr(l,1) = J_a; 
    
        J_a = J_a+dJ_a;

    end    
    
    I_0p50 = A*curr';
    T_0p50 = Temp_max';

    J_a = 0;

    J_a_max = 0.069;

    l = 0;    

    Dir = 'T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_0.75/Tx/';
    
    while J_a < J_a_max

        l = l+1;
        J_a_string = num2str(J_a);
    
        if length(J_a_string) == 7
            eval(['load ',Dir,J_a_string,'/ttx.dat']);
        elseif length(J_a_string) == 6
            eval(['load ',Dir,J_a_string,num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 5
            eval(['load ',Dir,J_a_string,num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 4
            eval(['load ',Dir,J_a_string,num2str(0),num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 3
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 2
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);            
        elseif length(J_a_string) == 1
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);            
        else
            eval(['load ',Dir,J_a_string,'/ttx.dat']);
        end
    
        Temp_max(l) = max(ttx(:,2));
        Temp_max(l) = T_0+Temp_max(l); % temperatura maxima no centro do vortice fechado
        curr(l,1) = J_a; 
    
        J_a = J_a+dJ_a;

    end    
    
    I_0p75 = A*curr';
    T_0p75 = Temp_max';    

    J_a = 0;

    J_a_max = 0.0705;

    l = 0;    

    Dir = 'T_0.75/kappa_1.0/Lz_2.4/hf_0.25_hs_1.00/Tx/';
    
    while J_a < J_a_max

        l = l+1;
        J_a_string = num2str(J_a);
    
        if length(J_a_string) == 7
            eval(['load ',Dir,J_a_string,'/ttx.dat']);
        elseif length(J_a_string) == 6
            eval(['load ',Dir,J_a_string,num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 5
            eval(['load ',Dir,J_a_string,num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 4
            eval(['load ',Dir,J_a_string,num2str(0),num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 3
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);
        elseif length(J_a_string) == 2
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);            
        elseif length(J_a_string) == 1
            eval(['load ',Dir,J_a_string,'.',num2str(0),num2str(0),num2str(0),num2str(0),num2str(0),'/ttx.dat']);            
        else
            eval(['load ',Dir,J_a_string,'/ttx.dat']);
        end
    
        Temp_max(l) = max(ttx(:,2));
        Temp_max(l) = T_0+Temp_max(l); % temperatura maxima no centro do vortice fechado
        curr(l,1) = J_a; 
    
        J_a = J_a+dJ_a;

    end    
    
    I_1p00 = A*curr';
    T_1p00 = Temp_max';

    %----------------------------------------------------------------------
    % \section{1.2 Faz o grafico comparativo entre as temperaturas
    %              no centro do vortice fechado}
    %----------------------------------------------------------------------
    
    figure;
    l_0p25 = plot(I_0p25,T_0p25,'LineWidth',1.5);
    hold on;
    l_0p50 = plot(I_0p50,T_0p50,'LineWidth',1.5);
    l_0p75 = plot(I_0p75,T_0p75,'LineWidth',1.5);
    l_1p00 = plot(I_1p00,T_1p00,'LineWidth',1.5);
    xlim([min([min(I_0p25) min(I_0p50) min(I_0p75) min(I_1p00)]) max([max(I_0p25) max(I_0p50) max(I_0p75) max(I_1p00)])]);
    ylim([min([min(T_0p25) min(T_0p50) min(T_0p75) min(T_1p00)]) max([max(T_0p25) max(T_0p50) max(T_0p75) max(T_1p00)])]);
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
    ylabel('$T/T_c$','Interpreter','LaTeX','FontSize',14);
    text(min(I_0p25)+0.1,0.95*max(T_0p25),'$T_0=0.75T_c,\,\kappa=1$','FontSize',24,'Interpreter','LaTeX');
    text(min(I_0p25)+0.1,0.935*max(T_0p25),'$l_z=2.4\xi(0)$','FontSize',24,'Interpreter','LaTeX');
    set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
    file_nt = strcat('T','_nt.png');
    file_t = strcat('T','_t.png');
    set(gca,'Color',[1 1 1]);
    background = get(gcf,'color');
    set(gcf,'color',[0.8 0.8 0.8]);
    set(gcf,'InvertHardCopy','off'); 
    print('-dpng', file_nt);
    cdata = imread(file_nt);
    imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
    
    cmd = ['mv -f *.png', ' ', 'T_0.75/kappa_1.0/Lz_2.4'];
    system(cmd);
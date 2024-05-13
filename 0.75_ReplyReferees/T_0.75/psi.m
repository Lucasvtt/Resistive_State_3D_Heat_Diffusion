close all;
clear all;

%--------------------------------------------------------------------------
% Dimensões da célula unitária
%--------------------------------------------------------------------------

ly = 16;
lz = 2.4;
dy = 0.4;
dz = 0.4;

y = 0:dy:ly;
z = 0:dz:lz;

%--------------------------------------------------------------------------
% Carrega os arquivos de psi
%--------------------------------------------------------------------------

load psiYZ.63;
psi_a = psiYZ;
load psiYZ.7;
psi_b = psiYZ;
load psiYZ.29;
psi_c = psiYZ;
load psiYZ.66;
psi_d = psiYZ;

%--------------------------------------------------------------------------
% Faz os gráficos de psi
%--------------------------------------------------------------------------

figure;
colormap('turbo');
psi_a = interp2(interp2(psi_a));
imagesc(y,z,psi_a,[0 max(max(psi_a))]);
axis equal;
axis tight;
axis xy;
xticks([]);
yticks([]);
set(gca,'visible','off')
% xlabel('$y/\xi$','Interpreter','LaTeX','FontSize',14);
% ylabel('$z/\xi$','Interpreter','LaTeX','FontSize',14);
set(gca,'TickLabelInterpreter','LaTeX','FontSize',17);
set(gca,'LineWidth',2);
set(gca,'TickDir','in');
%text(0.3,1.5,'(a)','Interpreter','LaTeX','FontSize',28,'Color','k');
set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
file_nt = strcat('psi_a','_nt.png');
file_t = strcat('psi_a','_t.png');
set(gca,'Color',[1 1 1]);
background = get(gcf,'color');
set(gcf,'color',[0.8 0.8 0.8]);
set(gcf,'InvertHardCopy','off'); 
print('-dpng', file_nt);
cdata = imread(file_nt);
imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);

figure;
colormap('turbo');
psi_b = interp2(interp2(psi_b));
imagesc(y,z,psi_b,[0 max(max(psi_b))]);
axis equal;
axis tight;
axis xy;
xticks([]);
yticks([]);
set(gca,'visible','off');
% xlabel('$y/\xi$','Interpreter','LaTeX','FontSize',14);
% ylabel('$z/\xi$','Interpreter','LaTeX','FontSize',14);
set(gca,'TickLabelInterpreter','LaTeX','FontSize',17);
set(gca,'LineWidth',2);
set(gca,'TickDir','in');
%text(0.3,1.5,'(b)','Interpreter','LaTeX','FontSize',28,'Color','k');
set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
file_nt = strcat('psi_b','_nt.png');
file_t = strcat('psi_b','_t.png');
set(gca,'Color',[1 1 1]);
background = get(gcf,'color');
set(gcf,'color',[0.8 0.8 0.8]);
set(gcf,'InvertHardCopy','off'); 
print('-dpng', file_nt);
cdata = imread(file_nt);
imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);

figure;
colormap('turbo');
psi_c = interp2(interp2(psi_c));
imagesc(y,z,psi_c,[0 max(max(psi_c))]);
axis equal;
axis tight;
axis xy;
xticks([]);
yticks([]);
set(gca,'visible','off')
% xlabel('$y/\xi$','Interpreter','LaTeX','FontSize',14);
% ylabel('$z/\xi$','Interpreter','LaTeX','FontSize',14);
set(gca,'TickLabelInterpreter','LaTeX','FontSize',17);
set(gca,'LineWidth',2);
set(gca,'TickDir','in');
%text(0.3,1.5,'(c)','Interpreter','LaTeX','FontSize',28,'Color','k');
set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
file_nt = strcat('psi_c','_nt.png');
file_t = strcat('psi_c','_t.png');
set(gca,'Color',[1 1 1]);
background = get(gcf,'color');
set(gcf,'color',[0.8 0.8 0.8]);
set(gcf,'InvertHardCopy','off'); 
print('-dpng', file_nt);
cdata = imread(file_nt);
imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);

figure;
colormap('turbo');
psi_d = interp2(interp2(psi_d));
imagesc(y,z,psi_d,[0 max(max(psi_d))]);
axis equal;
axis tight;
axis xy;
xticks([]);
yticks([]);
set(gca,'visible','off')
% xlabel('$y/\xi$','Interpreter','LaTeX','FontSize',14);
% ylabel('$z/\xi$','Interpreter','LaTeX','FontSize',14);
set(gca,'TickLabelInterpreter','LaTeX','FontSize',17);
set(gca,'LineWidth',2);
set(gca,'TickDir','in');
%text(0.3,1.5,'(d)','Interpreter','LaTeX','FontSize',28,'Color','k');
set(gca,'TickLabelInterpreter','LaTeX','FontSize',24);
file_nt = strcat('psi_d','_nt.png');
file_t = strcat('psi_d','_t.png');
set(gca,'Color',[1 1 1]);
background = get(gcf,'color');
set(gcf,'color',[0.8 0.8 0.8]);
set(gcf,'InvertHardCopy','off'); 
print('-dpng', file_nt);
cdata = imread(file_nt);
imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
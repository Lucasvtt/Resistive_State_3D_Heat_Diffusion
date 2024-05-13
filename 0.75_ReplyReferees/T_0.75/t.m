close all;
clear all;

%--------------------------------------------------------------------------
% Dimensões da célula unitária
%--------------------------------------------------------------------------

ly = 16;
lz = 2.4;
dy = 0.4;
dz = 0.4;
Ny = round(ly/dy);
Nz = round(lz/dz);

y = 0:dy:ly;
z = 0:dz:lz;

%--------------------------------------------------------------------------
% Carrega os arquivos de Temperatura no plano y
%--------------------------------------------------------------------------

load tyz.42;
tyz(1,1) = 0.5*(tyz(1,2)+tyz(2,1));
tyz(Nz+1,1) = 0.5*(tyz(Nz,1)+tyz(Nz+1,2));
tyz(1,Ny+1) = 0.5*(tyz(2,Ny+1)+tyz(1,Ny));
tyz(Nz+1,Ny+1) = 0.5*(tyz(Nz+1,Ny)+tyz(Nz,Ny+1));
tyz_a = tyz;

load tyz.43;
tyz(1,1) = 0.5*(tyz(1,2)+tyz(2,1));
tyz(Nz+1,1) = 0.5*(tyz(Nz,1)+tyz(Nz+1,2));
tyz(1,Ny+1) = 0.5*(tyz(2,Ny+1)+tyz(1,Ny));
tyz(Nz+1,Ny+1) = 0.5*(tyz(Nz+1,Ny)+tyz(Nz,Ny+1));
tyz_b = tyz;

load tyz.44;
tyz(1,1) = 0.5*(tyz(1,2)+tyz(2,1));
tyz(Nz+1,1) = 0.5*(tyz(Nz,1)+tyz(Nz+1,2));
tyz(1,Ny+1) = 0.5*(tyz(2,Ny+1)+tyz(1,Ny));
tyz(Nz+1,Ny+1) = 0.5*(tyz(Nz+1,Ny)+tyz(Nz,Ny+1));
tyz_c = tyz;

load tyz.45;
tyz(1,1) = 0.5*(tyz(1,2)+tyz(2,1));
tyz(Nz+1,1) = 0.5*(tyz(Nz,1)+tyz(Nz+1,2));
tyz(1,Ny+1) = 0.5*(tyz(2,Ny+1)+tyz(1,Ny));
tyz(Nz+1,Ny+1) = 0.5*(tyz(Nz+1,Ny)+tyz(Nz,Ny+1));
tyz_d = tyz;


%min_value = min([min(tyz_a(:)), min(tyz_b(:)), min(tyz_c(:)), min(tyz_d(:))]);
%max_value = max([max(tyz_a(:)), max(tyz_b(:)), max(tyz_c(:)), max(tyz_d(:))]);


tyz_max = max([max(max(tyz_a)) max(max(tyz_b)) max(max(tyz_c)) max(max(tyz_d))]);

%--------------------------------------------------------------------------
% Faz os gráficos de T
%--------------------------------------------------------------------------

figure;
colormap('hot');
tyz_a = interp2(interp2(tyz_a));
imagesc(y,z,tyz_a,[0 tyz_max]);
%imagesc(y, z, tyz_a, [min_value, max_value]);
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
file_nt = strcat('tyz_a','_nt.png');
file_t = strcat('tyz_a','_t.png');
set(gca,'Color',[1 1 1]);
background = get(gcf,'color');
set(gcf,'color',[0.8 0.8 0.8]);
set(gcf,'InvertHardCopy','off'); 
print('-dpng', file_nt);
cdata = imread(file_nt);
imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);

figure;
colormap('hot');
tyz_b = interp2(interp2(tyz_b));
imagesc(y,z,tyz_b,[0 tyz_max]);
%imagesc(y, z, tyz_b, [min_value, max_value]);
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
file_nt = strcat('tyz_b','_nt.png');
file_t = strcat('tyz_b','_t.png');
set(gca,'Color',[1 1 1]);
background = get(gcf,'color');
set(gcf,'color',[0.8 0.8 0.8]);
set(gcf,'InvertHardCopy','off'); 
print('-dpng', file_nt);
cdata = imread(file_nt);
imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);

figure;
colormap('hot');
tyz_c = interp2(interp2(tyz_c));
imagesc(y,z,tyz_c,[0 tyz_max]);
%imagesc(y, z, tyz_c, [min_value, max_value]);
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
file_nt = strcat('tyz_c','_nt.png');
file_t = strcat('tyz_c','_t.png');
set(gca,'Color',[1 1 1]);
background = get(gcf,'color');
set(gcf,'color',[0.8 0.8 0.8]);
set(gcf,'InvertHardCopy','off'); 
print('-dpng', file_nt);
cdata = imread(file_nt);
imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);

figure;
colormap('hot');
tyz_d = interp2(interp2(tyz_d));
imagesc(y,z,tyz_d,[0 tyz_max]);
%imagesc(y, z, tyz_d, [min_value, max_value]);
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
file_nt = strcat('tyz_d','_nt.png');
file_t = strcat('tyz_d','_t.png');
set(gca,'Color',[1 1 1]);
background = get(gcf,'color');
set(gcf,'color',[0.8 0.8 0.8]);
set(gcf,'InvertHardCopy','off'); 
print('-dpng', file_nt);
cdata = imread(file_nt);
imwrite(cdata,file_t,'png','BitDepth', 16,'transparency',[0.8 0.8 0.8]);
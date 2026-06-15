clearvars; close all; clc
set(groot,'DefaultFigurePosition', [100 100 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',10)
set(groot,'defaultaxesfontsize',18)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
list_factory = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_factory,'Interpreter'));for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
% MLF
addpath('/Users/charles/Documents/GIT/mlf')
%
VIEW    = [-160,40];
FSZ     = 16;
%%% Environment parameters
% Bounds
freq_bnd                = [1e-2 1e8];
porosity_bnd            = [.6 .99];
pore_mean_size_bnd      = [1e-6 1e-2];
pore_standard_dev_bnd   = [0 .5];

load('alpha.mat'), 
alpha_ref   = H;
alpha_app   = r;
load('beta.mat'),
beta_ref    = H;
beta_app    = r;

ZAref   = @(x) fun.impedence_abs(alpha_ref,beta_ref,x(:,1), x(:,2), x(:,3), x(:,4));
ZAapp   = @(x) fun.impedence_abs(alpha_app,beta_app,x(:,1), x(:,2), x(:,3), x(:,4));

%%% Plot some results
N   = [51 50 20 4];
x1  = 2*pi*logspace(log10(freq_bnd(1)),log10(freq_bnd(2)),N(1))*(1+rand(1)/50);
x2  = linspace(porosity_bnd(1),porosity_bnd(2),N(2))*(1+rand(1)/50);
x3  = logspace(log10(pore_mean_size_bnd(1)),log10(pore_mean_size_bnd(2)),N(3))*(1+rand(1)/50);
x4  = linspace(pore_standard_dev_bnd(1),pore_standard_dev_bnd(2),N(4))*(1+rand(1)/50);

%%%
x1label = '$x_1=\omega$ [rad/s]';
x2label = '$x_2=\phi$';
[X,Y]   = meshgrid(x1,x2);
h       = figure;
for i4 = 1:length(x4)
    for i3 = 1:length(x3)
        for jj = 1:numel(x2)
            for ii = 1:length(x1)
                pp              = [x3(i3) x4(i4)];
                p               = [1i*x1(ii) x2(jj) pp];
                tab_ref         = ZAref(p);
                tab_app         = ZAapp(p);
                %
                tab_ref_Z(jj,ii) = tab_ref(1);
                tab_ref_A(jj,ii) = tab_ref(2);
                tab_app_Z(jj,ii) = tab_app(1);
                tab_app_A(jj,ii) = tab_app(2);
            end
        end
        %
        clf,
        subplot(3,2,1); hold on, grid on
        surf(X,Y,tab_app_A,'EdgeColor','none'), hold on
        surf(X,Y,tab_ref_A,'EdgeColor','k','FaceColor','none')
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        zlabel('Absorption','Interpreter','latex')
        set(gca,'XScale','log'); 
        title('Absorbtion','Interpreter','latex')
        axis tight, view(VIEW(1),VIEW(2))
        %zlim([0 .6])
        subplot(3,2,2); hold on, grid on, axis tight
        imagesc(log10(abs(tab_ref_A-tab_app_A)/max(abs(tab_ref_A(:)))),'XData',x1,'YData',x2)
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        set(gca,'XScale','log');
        title('{\bf log}(abs. err.)/max.','Interpreter','latex')
        colorbar,
        %
        subplot(3,2,3); hold on, grid on
        surf(X,Y,real(tab_app_Z),'EdgeColor','none'), hold on
        surf(X,Y,real(tab_ref_Z),'EdgeColor','k','FaceColor','none')
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        zlabel('Impedance (real)','Interpreter','latex')
        set(gca,'XScale','log'); 
        title('Absorbtion','Interpreter','latex')
        axis tight, view(VIEW(1),VIEW(2))
        %
        subplot(3,2,4); hold on, grid on, axis tight
        imagesc(log10(abs(tab_ref_Z-tab_app_Z)/max(abs(tab_ref_A(:)))),'XData',x1,'YData',x2)
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        set(gca,'XScale','log');
        title('{\bf log}(abs. err.)/max.','Interpreter','latex')
        colorbar,
        %
        %
        subplot(3,2,5); hold on, grid on
        surf(X,Y,imag(tab_app_Z),'EdgeColor','none'), hold on
        surf(X,Y,imag(tab_ref_Z),'EdgeColor','k','FaceColor','none')
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        zlabel('Impedance (imag.)','Interpreter','latex')
        set(gca,'XScale','log'); 
        title('Absorbtion','Interpreter','latex')
        axis tight, view(VIEW(1),VIEW(2))
        %
        sgtitle({regexprep(func2str(H),'_','-'); ['$[x_3,x_4]=[\sigma_r,\overline{r}]=[' regexprep(num2str(pp),'\s*',',') ']$']},'interpreter','latex','FontSize',FSZ);
        drawnow
        %pause
        %fun.saveGIF(h,kk,name)
    end
end


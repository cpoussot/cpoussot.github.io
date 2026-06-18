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
CAS     = 1;
FSZ     = 16;
%%% Bounds
w_bnd                   = [1e-2 1e8];
porosity_bnd            = [.6 .99];
pore_mean_size_bnd      = [1e-6 1e-2];
pore_standard_dev_bnd   = [0 .5];
%%% Interpolation points
ip{1,1} = logspace(log10(w_bnd(1)),log10(w_bnd(2)),50);
ip{2,1} = linspace(porosity_bnd(1),porosity_bnd(2),10);
ip{3,1} = logspace(log10(pore_mean_size_bnd(1)),log10(pore_mean_size_bnd(2)),20);
ip{4,1} = linspace(pore_standard_dev_bnd(1),pore_standard_dev_bnd(2),20);
n = length(ip);
for ii = 1:n
    p_c{ii} = ip{ii}(2:2:end);
    p_r{ii} = ip{ii}(1:2:end);
end
%%% Define the alpha and beta functions
H_alpha = @(x) fun.dynamic_viscous_tortuosity(x(:,1), x(:,2), x(:,3), x(:,4));
H_beta  = @(x) fun.dynamic_thermal_compressibility(x(:,1), x(:,2), x(:,3), x(:,4));
%%% Data tensor
[y,x,dim]   = mlf.make_tab_vec(H_alpha,p_c,p_r);
tab_alpha   = mlf.vec2mat(y,dim);
[y,x,dim]   = mlf.make_tab_vec(H_beta,p_c,p_r);
tab_beta    = mlf.vec2mat(y,dim);

% N   = [1e4 3 1 1];
% x1  = logspace(log10(freq_bnd(1)),log10(freq_bnd(2)),N(1))*(1+rand(1)/50);
% x2  = 0.7;%linspace(porosity_bnd(1),porosity_bnd(2),N(2))*(1-rand(1)/50);
% x3  = [1e-4 5e-4 1e-3];%logspace(log10(pore_mean_size_bnd(1)),log10(pore_mean_size_bnd(2)),N(3))*(1+rand(1)/50);
% x4  = .1;%linspace(pore_standard_dev_bnd(1),pore_standard_dev_bnd(2),N(4))*(1+rand(1)/50);
% k = 0;
% figure, hold on
% for i4 = 1:length(x4)
%     for i3 = 1:length(x3)
%         for jj = 1:numel(x2)
%             for ii = 1:length(x1)
%                 tab_ref(ii,:) = ZAref([1i*x1(ii) x2(jj) x3(i3) x4(i4)]);
%             end
%             %plot(x1,tab_ref(:,2))
%             subplot(211), hold on
%             plot(x1,real(tab_ref(:,1)))
%             set(gca,'XScale','log')
%             %xlim([20 20000])
%             subplot(212), hold on
%             plot(x1,imag(tab_ref(:,1)))
%             set(gca,'XScale','log')%,'YScale','log')
%             %ylim([-1 1]*20)
%             %xlim([20 20000])
%             drawnow
%         end
%     end
% end
% %tab_app = ZAapp(p);

%%% Alg. 1: direct pLoe [A/G/P-V, 2025]
opt.method_null = 'svd0';
opt.method      = 'full';
opt.ord_obj     = [9 1 5 6];
opt.data_min    = false;
[r_alpha,imlf]  = mlf.alg1(tab_alpha,p_c,p_r,opt);
titre_alpha     = ['mLF alg. 1, $r=[' regexprep(num2str(imlf.ord),'\s*',',') ']$']
opt.ord_obj     = [8 1 6 5];
[r_beta,imlf]   = mlf.alg1(tab_beta,p_c,p_r,opt);
titre_beta      = ['mLF alg. 1, $r=[' regexprep(num2str(imlf.ord),'\s*',',') ']$']

%%
H = H_alpha; r = r_alpha; titre = titre_alpha; name = 'alpha';
%H = H_beta;  r = r_beta; titre = titre_beta; name = 'beta';

%%% Plot some results
N   = [51 50 20 3];
x1  = logspace(log10(w_bnd(1)),log10(w_bnd(2)),N(1))*(1+rand(1)/50);
x2  = linspace(porosity_bnd(1),porosity_bnd(2),N(2))*(1+rand(1)/50);
x3  = logspace(log10(pore_mean_size_bnd(1)),log10(pore_mean_size_bnd(2)),N(3))*(1+rand(1)/50);
x4  = linspace(pore_standard_dev_bnd(1),pore_standard_dev_bnd(2),N(4))*(1+rand(1)/50);

%%%
x1label = '$x_1=\omega$ [rad/s]';
x2label = '$x_2=\sigma_r$';
[X,Y]   = meshgrid(x1,x2);
h       = figure('Color','white');
kk      = 0;
for i4 = 1:length(x4)
    for i3 = 1:length(x3)
        for jj = 1:numel(x2)
            for ii = 1:length(x1)
                pp              = [x3(i3) x4(i4)];
                p               = [1i*x1(ii) x2(jj) pp];
                tab_ref         = H(p);
                tab_app1        = r(p);
                %
                tab_refR(jj,ii) = real(tab_ref);
                tab_refI(jj,ii) = imag(tab_ref);
                tab_app1R(jj,ii)= real(tab_app1);
                tab_app1I(jj,ii)= imag(tab_app1);
            end
        end
        %
        clf,
        subplot(2,2,1); hold on, grid on
        surf(X,Y,tab_app1R,'EdgeColor','none'), hold on
        surf(X,Y,tab_refR,'EdgeColor','k','FaceColor','none')
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        zlabel('Real(.)','Interpreter','latex')
        set(gca,'XScale','log'); 
        legend({'Rational approximation' 'Original model'},'Location','NorthWest')
        title(titre,'Interpreter','latex')
        axis tight, view(VIEW(1),VIEW(2))
        %zlim([min(tab_refR(:)) max(tab_refR(:))])
        zlim([0 .6])
        subplot(2,2,2); hold on, grid on, axis tight
        imagesc(log10(abs(tab_refR-tab_app1R)/max(abs(tab_refR(:)))),'XData',x1,'YData',x2)
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        set(gca,'XScale','log');
        title('{\bf log}(abs. err./max.)','Interpreter','latex')
        colorbar,
        %clim([-12 0])
        %
        subplot(2,2,3); hold on, grid on
        surf(X,Y,tab_app1I,'EdgeColor','none'), hold on
        surf(X,Y,tab_refI,'EdgeColor','k','FaceColor','none')
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        zlabel('Imag(.)','Interpreter','latex')
        set(gca,'XScale','log'); 
        title(titre,'Interpreter','latex')
        axis tight, view(VIEW(1),VIEW(2))
        %zlim([min(tab_refI(:)) max(tab_refI(:))]), 
        zlim([-.175 0])
        subplot(2,2,4); hold on, grid on, axis tight
        imagesc(log10(abs(tab_refI-tab_app1I)/max(abs(tab_refI(:)))),'XData',x1,'YData',x2)
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        set(gca,'XScale','log');
        title('{\bf log}(abs. err./max.)','Interpreter','latex')
        colorbar,
        %clim([-12 0])
        sgtitle({regexprep(func2str(H),'_','-'); ['$[x_3,x_4]=[\phi,\overline{r}]=[' regexprep(num2str(pp),'\s*',',') ']$']},'interpreter','latex','FontSize',FSZ);
        drawnow
        kk = kk + 1; fun.saveGIF(h,kk,name)
    end
end
%%
ZAref   = @(x) fun.impedence_abs(H_alpha,H_beta, x(:,1), x(:,2), x(:,3), x(:,4));
ZAapp   = @(x) fun.impedence_abs(r_alpha,r_beta, x(:,1), x(:,2), x(:,3), x(:,4));

x1      = logspace(log10(w_bnd(1)),log10(w_bnd(2)),1e4)*(1+rand(1)/50);
x2      = 0.7;
x3      = logspace(-6,-3,40)*(1+rand(1)/100);
x4      = .1;
kk      = 0;
col     = hsv(5);
h=figure('Color','white'), hold on, grid on, axis tight
for i4 = 1:length(x4)
    for i3 = 1:length(x3)
        for i2 = 1:numel(x2)
            for i1 = 1:length(x1)
                Abs_ref(i1,:) = ZAref([1i*x1(i1) x2(i2) x3(i3) x4(i4)]);
                Abs_app(i1,:) = ZAapp([1i*x1(i1) x2(i2) x3(i3) x4(i4)]);
            end
            cla
            h1=plot(x1,Abs_ref(:,2),'LineWidth',3);%,'Color',col(kk,:));
            h2=plot(x1,Abs_app(:,2),'k--','LineWidth',3);
            set(gca,'XScale','log')
            xlabel('Frequency [rad/s]')
            ylabel('Absorption coefficient')
            ylim([0 1])
            title(['$\{\sigma_r,\phi,\overline{r}\}=\{' num2str(x2(i2),2) ',' num2str(x3(i3),2) ',' num2str(x4(i4),2) '\}$'])
            legend({'Original model' 'Rational approximation'},'Location','west')
            drawnow
            kk = kk + 1; fun.saveGIF(h,kk,'absorption',.2)
        end
    end
end

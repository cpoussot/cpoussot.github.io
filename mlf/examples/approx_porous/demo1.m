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
%%% Environment parameters
% Bounds
freq_bnd                = [1e-2 1e8];
porosity_bnd            = [.6 .99];
pore_mean_size_bnd      = [1e-6 1e-2];
pore_standard_dev_bnd   = [0 .5];
% Material parameters variables
% omega               = 2*pi*logspace(log10(freq_bnd(1)),log10(freq_bnd(2)),50);
% porosity            = linspace(porosity_bnd(1),porosity_bnd(2),10);
% pore_mean_size      = linspace(pore_mean_size_bnd(1),pore_mean_size_bnd(2),20);
% pore_standard_dev   = linspace(pore_standard_dev_bnd(1),pore_standard_dev_bnd(2),20);

omega               = 2*pi*logspace(log10(freq_bnd(1)),log10(freq_bnd(2)),50);
porosity            = linspace(porosity_bnd(1),porosity_bnd(2),10);
pore_mean_size      = logspace(log10(pore_mean_size_bnd(1)),log10(pore_mean_size_bnd(2)),20);
pore_standard_dev   = linspace(pore_standard_dev_bnd(1),pore_standard_dev_bnd(2),20);

%%% Functions, ip
ip{1,1} = omega;
ip{2,1} = porosity;
ip{3,1} = pore_mean_size;
ip{4,1} = pore_standard_dev;
ii      = 1;
% % iR => conj()
% i1      = 2;
% p_c{1}  = -1i*omega(2:2:end);
% p_r{1}  = -1i*omega(1:2:end);
% p_c{1}  = sort([p_c{1} conj(p_c{1})]);
% p_r{1}  = sort([p_r{1} conj(p_r{1})]);
% ip{1}   = [p_c{1} p_r{1}];
% R
%%% Data tensor
n = length(ip);
for ii = ii:n
    p_c{ii} = ip{ii}(2:2:end);
    p_r{ii} = ip{ii}(1:2:end);
end
[y,x,dim]   = mlf.make_tab_vec(H,p_c,p_r);
tab         = mlf.vec2mat(y,dim);

switch CAS
    case 1
        H       = @(x) fun.dynamic_viscous_tortuosity(x(:,1), x(:,2), x(:,3), x(:,4));
        name    = 'alpha';
        %
        METH    = 'full';
        ord_tol = 1e-9;
        %ord_obj = [18 1 5 6];
        ord_obj = [9 1 5 6];
        %ord_obj = [18 inf inf inf];
    case 2
        H       = @(x) fun.dynamic_thermal_compressibility(x(:,1), x(:,2), x(:,3), x(:,4));
        name    = 'beta';
        %
        METH    = 'full';
        ord_tol = 1e-10;
        ord_obj = [13 1 6 5];
        ord_obj = [8 1 6 5];
    % case 3
    %     H       = @(x) fun.impedence_abs(x(:,1), x(:,2), x(:,3), x(:,4));
    %     METH    = 'full';
    %     ord_tol = 1e-5;
    %     %ord_tol = 1e-9;
    %     ord_obj = [13 1 6 5];
end

%%% Alg. 1: direct pLoe [A/G/P-V, 2025]
opt = [];
tic;
opt.method_null = 'svd0';
opt.method      = METH;
opt.ord_tol     = ord_tol;
opt.ord_obj     = ord_obj;
opt.ord_N       = 10;
opt.ord_show    = true;
opt.data_min    = false;
[r,imlf]        = mlf.alg1(tab,p_c,p_r,opt);
if opt.ord_show; drawnow, mlf.figSavePDF(['svd_' num2str(CAS)],.5), pause(.5), end
toc
titre = ['mLF alg. 1, $r=[' regexprep(num2str(imlf.ord),'\s*',',') ']$'];
save(name,'H','r')

% %%% Realization Lagrangian
% % Original
% [Glag,ilagREAL] = mlf.make_realization_lag(imlf.pc,imlf.w,imlf.c,[]);
% [Glag,ilagREAL] = mlf.make_realization_compressed(ilagREAL);
% %%% Monomial real
% opt_real        = [];
% opt_real.s_gam  = 1;
% opt_real.s_del  = 2:n;
% [Gmon,imonREAL] = mlf.make_realization_compressed_mon(imlf.pc,imlf.w,imlf.c,opt_real);
%%
%r  = @(x) Glag(x(:,1),x(:,2),x(:,3),x(:,4));
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
        title(titre,'Interpreter','latex')
        axis tight, view(VIEW(1),VIEW(2))
        %zlim([min(tab_refR(:)) max(tab_refR(:))])
        zlim([0 .6])
        subplot(2,2,2); hold on, grid on, axis tight
        imagesc(log10(abs(tab_refR-tab_app1R)/max(abs(tab_refR(:)))),'XData',x1,'YData',x2)
        xlabel(x1label,'Interpreter','latex')
        ylabel(x2label,'Interpreter','latex')
        set(gca,'XScale','log');
        title('{\bf log}(abs. err.)/max.','Interpreter','latex')
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
        title('{\bf log}(abs. err.)/max.','Interpreter','latex')
        colorbar,
        %clim([-12 0])
        sgtitle({regexprep(func2str(H),'_','-'); ['$[x_3,x_4]=[\sigma_r,\overline{r}]=[' regexprep(num2str(pp),'\s*',',') ']$']},'interpreter','latex','FontSize',FSZ);
        drawnow
        %fun.saveGIF(h,kk,name)
    end
end


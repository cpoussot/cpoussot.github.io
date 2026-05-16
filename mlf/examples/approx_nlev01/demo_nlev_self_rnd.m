clearvars; close all; format compact; format long e; clc;
set(groot,'DefaultFigurePosition', [300 300 1000 400]);
set(groot,'defaultlinelinewidth',5)
set(groot,'defaultlinemarkersize',20)
set(groot,'defaultaxesfontsize',24)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
% mLF
addpath('/Users/charles/Documents/GIT/mLF')
addpath('/Users/charles/Documents/GIT/LF')
%
SAVEIT  = true;
%%% Examples article
E   = diag(logspace(-4,10,10));
Phi = @(x) (x(:,1)+0.01*exp(-x(:,1).*x(:,2)))*eye(10)+E;
H   = @(x) ones(1,10)*(Phi(x)\ones(10,1) );
Hf  = @(x1,x2) H([x1,x2]);
% IP
rng(1)
ip{1} = .2*(rand(22,1)-.5+1i*.5*rand(22,1)); % then complex conjugated
ip{1} = ip{1}(2:end-1);
ip{2} = 18:52;%
% interlace
for ii = 1:numel(ip)
    p_c{ii} = ip{ii}(2:2:end);
    p_r{ii} = ip{ii}(1:2:end);
end
% complex conjugate ip{1}
pc{1} = [];
pr{1} = [];
for ii = 1:length(p_c{1})
    pc{1} = [pc{1} p_c{1}(ii) conj(p_c{1}(ii))];
    pr{1} = [pr{1} p_r{1}(ii) conj(p_r{1}(ii))];
end
p_c{1}  = [pc{1}];
p_r{1}  = [pr{1}];
ip{1}   = [p_c{1} p_r{1}];
%%% Tensor
tab = mlf.make_tab(Hf,p_c,p_r,true);
%%% Lagrangian multivariate
opt.ord_tol     = 1e-12;
opt.method_null = 'svd';
opt.method      = 'full';
opt.ord_obj     = [inf 8];
opt.ord_show    = true;
[g,iloe]        = mlf.alg1(tab,p_c,p_r,opt);
if SAVEIT; drawnow, mlf.figSavePNG('svd_rnd',.5), pause(.5); end
%%% Realization Lagrangian
% Original
[~,ireal]   = mlf.make_realization_lag(iloe.pc,iloe.w,iloe.c,[]);
% Compressed
[H,ireal]   = mlf.make_realization_compressed(ireal);
%%
C   = .075*exp(1i*linspace(0,2*pi,1e3));
handler = figure; hold on, grid on
set(gcf, 'Color', 'white')
pSpace  = ip{2};%linspace(min(ip{2}),max(ip{2}),50);
pSpace  = linspace(min(ip{2}),max(ip{2}),80);
%pSpace  = [20 30 35 50];
for ii = 1:numel(pSpace)
    cla
    p   = pSpace(ii);
    Hfp = @(x) Hf(x,p);    
    %%% Loewner classic
    la  = p_c{1}; k = length(la); R = ones(1,k);
    mu  = p_r{1}; q = length(mu); L = ones(q,1);
    for i = 1:k; W(1,1,i) = Hfp(la(i)); end
    for i = 1:q; V(1,1,i) = Hfp(mu(i)); end
    [hloe,info_loe] = lf.loewner_tng(la,mu,W,V,R,L);
    lam_classik = eig(info_loe.Hr);
    % Phi = sE-A
    Phir        = ireal.Phi;
    A           = -Phir(0,p);
    E           = Phir(1,p)+A;
    [eigV,eigv] = eig(A,E); eigv = diag(eigv); 
    eigv(isinf(eigv))=[];
    eigv(isnan(eigv))=[];
    uns = numel(find(eigv(real(eigv)>0)));
    %
    %plot(real(C),imag(C),'k--','DisplayName','Contour $\partial\Omega$')
    plot(real(ip{1}),imag(ip{1}),'.','DisplayName','$z(1,\cdots,n_1)$') 
    plot(real(iloe.pc{1}),imag(iloe.pc{1}),'s','DisplayName','$\lambda_1$') 
    plot(real(iloe.pr{1}),imag(iloe.pr{1}),'d','DisplayName','$\mu_1$') 
    %plot(real(lam_classik),imag(lam_classik),'^','DisplayName','$\lambda$ (std. Loewner, $p$ frozen)')
    plot(real(eigv),imag(eigv),'v','DisplayName','$\lambda$ (est.)')
    xlabel('$\textrm{Re}(z)$','Interpreter','latex')
    ylabel('$\textrm{Im}(z)$','Interpreter','latex')
    title(sprintf('$p=%.2f$ (%d unstable)',p,uns))
    legend('show','Location','eastoutside')
    set(gca,'xlim',.15*[-1 1],'ylim',.15*[-1 1])
    axis square
    drawnow, %pause
    %if SAVEIT; mlf.figSavePDF(['figures/ex_' num2str(2) '_' num2str(ii)]); end
    if SAVEIT; mlf.saveGIF(handler,ii,['ex_rnd_' num2str(2) '_min' num2str(min(pSpace)) 'max_' num2str(max(pSpace)) ]); end
end
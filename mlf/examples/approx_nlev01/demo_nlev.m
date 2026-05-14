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

%%% Examples article
SAVEIT  = true;
%EX_NUM  = 1; tol = 1e-8; ord = [];
EX_NUM  = 2; tol = 1e-12; ord = [inf 8];
%EX_NUM  = 3; tol = 1e-5;
[Hf,Phi,ip,p_c,p_r,C,S,XLIM,YLIM]   = examples(EX_NUM);
tab                                 = mlf.make_tab(Hf,p_c,p_r,true);

%%% Lagrangian multivariate
opt.ord_tol     = tol;
opt.method_null = 'svd0';
opt.method      = 'full';
opt.ord_obj     = ord;
opt.ord_show    = true;
opt.data_min    = true;
[Gloe,iloe]     = mlf.alg1(tab,p_c,p_r,opt);
if SAVEIT; drawnow, mlf.figSavePNG('svd',.5), pause(.5); end

%%% Realization Lagrangian
% Original
[~,info_l]      = mlf.make_realization_lag(iloe.pc,iloe.w,iloe.c,[]);
% Compressed
[Hrl,info_l]    = mlf.make_realization_compressed(info_l);

%%
handler=figure, hold on, grid on
pSpace  = linspace(min(ip{2}),max(ip{2}),50);
%pSpace  = [20 30 35 50]
for ii = 1:numel(pSpace)
    %subplot(121)
    cla
    p   = pSpace(ii);
    Hfp = @(x) Hf(x,p);
    %%% Lambert
    
    %%% Loewner classic
    la  = p_c{1}; k = length(la); R = ones(1,k);
    mu  = p_r{1}; q = length(mu); L = ones(q,1);
    for i = 1:k; W(1,1,i) = Hfp(la(i)); end
    for i = 1:q; V(1,1,i) = Hfp(mu(i)); end
    %opt.target = 40;
    [hloe,info_loe] = lf.loewner_tng(la,mu,W,V,R,L);
    lam_classik = eig(info_loe.Hr);
    % Phi = sE-A
    Phir        = info_l.Phi;
    A           = -Phir(0,p);
    E           = Phir(1,p)+A;
    [eigV,eigv] = eig(A,E); eigv = diag(eigv);
    % for i = 1:length(lam)
    %     norm(V(:,i)*Phi([lam(i),p]) - V(:,i)*lam(i))
    % end
    % % TF
    % Htf         = mlf.tfp(iloe.pc,iloe.w,iloe.c,p);
    % Htf         = tf(minreal(Htf));
    % [b,a]       = tfdata(Htf);
    % [res,lam,k] = residue(b{1},a{1});
    %
    %plot(real(C),imag(C),'k--','DisplayName','Contour $\partial\Omega$')
    plot(real(ip{1}),imag(ip{1}),'.','DisplayName','$z_1(1,\cdots,n_1)$') 
    plot(real(iloe.pc{1}),imag(iloe.pc{1}),'s','DisplayName','$\lambda_1$') 
    plot(real(iloe.pr{1}),imag(iloe.pr{1}),'d','DisplayName','$\mu_1$') 
    plot(real(lam_classik),imag(lam_classik),'^','DisplayName','$\lambda$ (std. Loewner, $p$ frozen)')
    plot(real(eigv),imag(eigv),'v','DisplayName','$\lambda$ (est.)')
    plot(real(S(p)),imag(S(p)),'o','DisplayName','$\lambda$ (exact)')
    xlabel('$\textrm{Re}(z)$','Interpreter','latex')
    ylabel('$\textrm{Im}(z)$','Interpreter','latex')
    title(sprintf('$p=%.2f$',p))
    legend('show','Location','eastoutside')
    set(gca,'xlim',XLIM,'ylim',YLIM)
    axis square
    drawnow, %pause
    %if SAVEIT; mlf.figSavePDF(['figures/ex_' num2str(EX_NUM) '_' num2str(ii)]); end
    if SAVEIT; mlf.saveGIF(handler,ii,['ex_' num2str(EX_NUM) '_min' num2str(min(pSpace)) 'max_' num2str(max(pSpace)) ]); end
end
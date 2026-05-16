clearvars, close all, clc
set(groot,'DefaultFigurePosition', [300 100 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for iloe1 = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(iloe1)},'factory','default'),'latex'); end
addpath('/Users/charles/Documents/GIT/mlf')

%%% Define model
rng(1)
n       = 3;
H       = @(x) x(:,1).*x(:,2)+x(:,1).*x(:,3)+x(:,2).*x(:,3);

%%% Define interpolation points along each variables
ip{1}   = linspace(-1,1,10);
ip{2}   = linspace(-1,1,10);
ip{3}   = linspace(-1,1,10);

%%% Split interpolation points into column and row data sets and build tensor
for i = 1:n
    p_c{i} = ip{i}(2:2:end); 
    p_r{i} = ip{i}(1:2:end);
end
[y,x,dim]   = mlf.make_tab_vec(H,p_c,p_r);
T           = mlf.vec2mat(y,dim);

%%% Alg. 1: direct pLoe [A/G/P-V, 2025]
opt.ord_tol     = 1e-8;
opt.ord_show    = true;
[g,iloe]        = mlf.alg1(T,p_c,p_r,opt);
%drawnow, mlf.figSavePNG('svd',.5), pause(.5)

%%% Some mismatch comparison
for i = 1:1e3
    x_try   = rand(1,3);
    err(i)  = abs(H(x_try)-g(x_try));
end
mean(err)

%%% KST
[Bary,Lag,Cx]   = mlf.decoupling(iloe);
PHI1            = Bary{1}.*Lag{1};
PHI2            = Bary{2}.*Lag{2};
PHI3            = Bary{3}.*Lag{3};
num             = simplify(sum(iloe.w.*PHI1.*PHI2.*PHI3));
den             = simplify(sum(PHI1.*PHI2.*PHI3));
simplify(num/den)
phi1 = latex((PHI1));
phi2 = latex((PHI2));
phi3 = latex((PHI3));
for ii = 1:numel(p_c)
    phi1 = strrep(phi1,['s_{' num2str(ii) '}'],['x_{' num2str(ii) '}']);
    phi2 = strrep(phi2,['s_{' num2str(ii) '}'],['x_{' num2str(ii) '}']);
    phi3 = strrep(phi3,['s_{' num2str(ii) '}'],['x_{' num2str(ii) '}']);
end

%%% Plot some results
% Along first and second variables 
% Other variables are randomly chosen between bounds
x1      = linspace(min(ip{1}),max(ip{1}),40)+rand(1)/10;
x2      = linspace(min(ip{2}),max(ip{2}),41)+rand(1)/10;
[X,Y]   = meshgrid(x1,x2);
rnd_p   = mlf.rand(n-2,p_r(3:end),false);
for ii = 1:numel(x1)
    for jj = 1:numel(x2)
        param           = [x1(ii) x2(jj) rnd_p];
        tab_ref(jj,ii)  = H(param);
        tab_app1(jj,ii) = g(param);
    end
end
%
titre1 = ['mLF alg. 1 (direct), $r=[' regexprep(num2str(iloe.ord),'\s*',',') ']$'];
figure
subplot(1,2,1); hold on, grid on
surf(X,Y,tab_app1,'EdgeColor','none'), hold on
surf(X,Y,tab_ref,'EdgeColor','k','FaceColor','none')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title(titre1,'Interpreter','latex')
legend('Approximation $g$','Original $H$')
axis tight, zlim([min(tab_ref(:)) max(tab_ref(:))]), view(-20,15)
subplot(1,2,2); hold on, grid on, axis tight
imagesc(log10(abs(tab_ref-tab_app1)/max(abs(tab_ref(:)))),'XData',x1,'YData',x2)
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('{\bf log}(abs. err.)/max.','Interpreter','latex')
colorbar,
sgtitle('$x_1x_2+x_1x_3+x_2x_3$','Interpreter','latex','FontSize',25)
drawnow, mlf.figSavePNG('eval',.5)


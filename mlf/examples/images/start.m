clearvars, close all, clc
set(groot,'DefaultFigurePosition', [300 100 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for iloe1 = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(iloe1)},'factory','default'),'latex'); end
addpath('/Users/charles/Documents/GIT/mlf')

%%% Define model
n       = 2;
H       = @(x) sqrt(x(:,1)).*x(:,2)+sin(x(:,1)+x(:,2).^2);

%%% Define interpolation points along each variables
ip{1}   = linspace(0,2,20);
ip{2}   = linspace(-pi,pi,51);

%%% Split interpolation points into column and row data sets and build tensor
for i = 1:n
    p_c{i} = ip{i}(2:2:end); 
    p_r{i} = ip{i}(1:2:end);
end
[y,x,dim]   = mlf.make_tab_vec(H,p_c,p_r);
T           = mlf.vec2mat(y,dim);

%%% Alg. 1: direct pLoe [A/G/P-V, 2025]
opt.ord_tol = 1e-6;
[g,iloe]    = mlf.alg1(T,p_c,p_r,opt);

%%% Some mismatch comparison
for i = 1:1e2
    x_try   = rand(1,2);
    err(i)  = abs(H(x_try)-g(x_try));
end
mean(err)
 
%%% Plot some results
% Along first and second variables 
% Other variables are randomly chosen between bounds
x1      = linspace(min(ip{1}),max(ip{1}),40)+rand(1)/10;
x2      = linspace(min(ip{2}),max(ip{2}),41)+rand(1)/10;
[X,Y]   = meshgrid(x1,x2);
[Xd,Yd] = meshgrid(ip{1},ip{2});
rnd_p   = mlf.rand(n-2,p_r(3:end),false);
for ii = 1:numel(x1)
    for jj = 1:numel(x2)
        param           = [x1(ii) x2(jj) rnd_p];
        tab_ref(jj,ii)  = H(param);
        tab_app1(jj,ii) = g(param);
    end
end
%
orders = ['complexity along $(x_1,x_2)$ is $(' regexprep(num2str(iloe.ord),'\s*',',') ')$'];
figure
subplot(1,2,1); hold on, grid on
plot3(x(:,1),x(:,2),y,'k.','Color',[1 1 1]*.2), hold on
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('Data tensor','Interpreter','latex')
axis tight, zlim([min(tab_ref(:)) max(tab_ref(:))]), view(-20,15)
subplot(1,2,2); hold on, grid on, axis tight
plot3(x(:,1),x(:,2),y,'k.','Color',[1 1 1]*.2), hold on
surf(X,Y,tab_app1,'EdgeColor','none')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title({'Data tensor vs. Approximation $g(x_1,x_2)$'; orders},'Interpreter','latex')
axis tight, zlim([min(tab_ref(:)) max(tab_ref(:))]), view(-20,15)
%
sgtitle(['$' num2str(numel(y)) '$ tensorized data'],'Interpreter','latex','FontSize',25)
drawnow, mlf.figSavePNG('approx',.5)


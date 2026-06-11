clearvars, close all, clc
%%% Chose model
syms s1 s2
dt1     = .01;
dt2     = dt1/2;%2;
Hsym    = s2/(exp(dt1*s1)^2-1/2) + ...
           1/(exp(dt2*s1)+1/3);
H       = matlabFunction(Hsym);
H       = @(x) H(x(:,1),x(:,2));
%%% Define IP
ip{1}   = linspace(1e-3,pi/dt1,20);
ip{2}   = linspace(1e-2,1,10);
n       = numel(ip);
%%% Split IP and build tensor
for i = 1:n
    p_c{i} = ip{i}(2:2:end); 
    p_r{i} = ip{i}(1:2:end);
end
%%% Data tensor/rand
[y,x,dim]   = mlf.make_tab_vec(H,p_c,p_r);
tab         = mlf.vec2mat(y,dim);
%%% Alg. 1: direct pLoe [A/G/P-V, 2025]
opt.ord_tol     = 1e-8;  % SVD tolerance
opt.method_null = 'svd'; % null space method
opt.method      = 'rec'; % full or recursive method
opt.ord_show    = false;  % show order detection step
[glag,imlf]     = mlf.alg1(tab,p_c,p_r,opt);
%%% Monomial real
opt_real.s_gam  = 1;
opt_real.s_del  = 2;
[gmon,imon_tf]  = mlf.tf_monomial(imlf.pc,imlf.w,imlf.c,false);
[~,imon_re0]    = mlf.make_realization_mon(imlf.pc,imlf.w,imlf.c,opt_real);
[~,imon_re1]    = mlf.make_realization_compressed(imon_re0);
[Gmon,imon_red] = mlf.make_realization_compressed_mon(imlf.pc,imlf.w,imlf.c,opt_real);
%%% DSS
Phir0p          = subs(imon_red.Phi_s,'s1',0);
Phir1p          = subs(imon_red.Phi_s,'s1',1);
Er              = double(Phir1p-Phir0p);
%Ar              = matlabFunction(-Phir0p);
Ar              = double(-Phir0p);
Br              = double(imon_red.Gr);
Cr              = @(p) double(imon_red.Wr(0,p));%matlabFunction()
Gss             = @(p) dss(Ar,Br,Cr(p),0,Er);
vpa(imon_red.Phi_s,2)
vpa(Phir0p,2)
vpa(Phir1p,2)
% %%% KST
% [Bary,Lag,Cx]   = mlf.decoupling(imlf);
% PHI1            = Bary{1}.*Lag{1};
% PHI2            = Bary{2}.*Lag{2};
% num             = simplify(sum(imlf.w.*PHI1.*PHI2));
% den             = simplify(sum(PHI1.*PHI2));
% Hkst            = simplify(num/den);
% vpa(Hkst,3)
%
Hz = @(p) tf(p,[1 0 -1/2],dt1)+tf(1,[1 1/3],dt1) 
pSpan = linspace(min(ip{2}),max(ip{2}),10)+rand(1)/10;
tSpan_s = 0:dt1:.1;
for i = 1:numel(pSpan)
    y0(:,i) = step(Hz(pSpan(i)),tSpan_s);
end
tSpan_c = linspace(0,tSpan_s(end),1e3);
for i = 1:numel(pSpan)
    y1(:,i) = step(Gss(pSpan(i)),tSpan_c);
end

figure, hold on
plot(tSpan_s,y0)
plot(tSpan_c,y1,'k--')
% %%
% %%% rand
% for i = 1:1e3
%     x_try   = rand(1,n);
%     err(i)  = abs(H(x_try)-g(x_try));
% end
% mean(err)
% %%% Along first and second variables 
% x1      = logspace(log10(min(ip{1})),log10(max(ip{1})),60)+rand(1)/10;
% x2      = linspace(min(ip{2}),max(ip{2}),41)+rand(1)/10;
% [X,Y]   = meshgrid(x1,x2);
% rnd_p   = [];
% if n > 2; rnd_p = mlf.rand(n-2,p_r(3:end),false); end
% for ii = 1:numel(x1)
%     for jj = 1:numel(x2)
%         param           = [x1(ii) x2(jj) rnd_p];
%         paramStr        = regexprep(num2str(param,36),'\s*',',');
%         tab_ref(jj,ii)  = H(param);
%         tab_app(jj,ii)  = g(param);
%     end
% end
% %
% figure
% subplot(1,2,1); hold on, grid on
% surf(X,Y,abs(tab_app),'EdgeColor','none'), hold on
% surf(X,Y,abs(tab_ref),'EdgeColor','k','FaceColor','none')
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% title('Original vs. Approximation','Interpreter','latex')
% legend('Approximation $g$','Original $H$')
% axis tight,% zlim([min(tab_ref(:)) max(tab_ref(:))]), 
% view(-30,40)
% set(gca,'XScale','log')
% subplot(1,2,2); hold on, grid on, axis tight
% imagesc(log10(abs(tab_ref-tab_app)/max(abs(tab_ref(:)))),'XData',x1,'YData',x2)
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
% title('{\bf log}(abs. err.)/max.','Interpreter','latex')
% colorbar,
% %name = infoCas.name;
% %for ii = 1:numel(p_c), name = strrep(name,['\var{' num2str(ii) '}'],['x_{' num2str(ii) '}']); end
% %sgtitle(name,'Interpreter','latex','FontSize',25)
% %drawnow, mlf.figSavePNG('eval',.5), pause(.5)
% 
% 
% 
% 
% imon_red
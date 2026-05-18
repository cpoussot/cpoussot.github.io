clearvars; close all; format compact; format long e; clc;
set(groot,'DefaultFigurePosition', [300 300 1000 400]);
set(groot,'defaultlinelinewidth',5)
set(groot,'defaultlinemarkersize',20)
set(groot,'defaultaxesfontsize',24)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
%%% model
tau = .4;
CAS = 1
switch CAS
    case 1
        A0  = -2;
        A1  = -1;
        E   = eye(1);
        C   = 1;
        B   = C';
        H   = @(s) C*((E*s-A0-A1*exp(-tau*s))\B);
    case 2
        A0  = diag([1 10]);
        A1  = diag([-1e-2 -30]);
        E   = eye(length(A1));
        B   = ones(length(A1),1);
        C   = B';
        H   = @(s) C*((E*s-A0-A1*exp(-tau*s))\B);
end
%%% LF
% IP
Radius = 20;
ip{1} = Radius*exp(1i*linspace(0,pi,100)); % then complex conjugated
ip{1} = ip{1}(2:end-1);
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
%
la  = p_c{1}; k = length(la); R = ones(1,k);
mu  = p_r{1}; q = length(mu); L = ones(q,1);
for i = 1:k; W(1,1,i) = H(la(i)); end
for i = 1:q; V(1,1,i) = H(mu(i)); end
[hloe,info_loe] = lf.loewner_tng(la,mu,W,V,R,L);
lam_lf = eig(info_loe.Hr);

lam_lambert = lambert(A0,A1,tau); 
% lamAE   = eig(A,E);
% a1      = lamAE(:);
% a0      = diag(A0);
% for k = 0:10
%     for i = 1:numel(lamAE)
%         lam(i,k+1) = lambertw(k,tau*a1(i)*exp(-a0(i)*tau))/tau+a0(i);
%     end
% end
% lam_lambert = [lam conj(lam)].';
lam_lambert = lam_lambert(:);

figure, grid on, hold on
plot(real(ip{1}),imag(ip{1}),'.','DisplayName','$z(1,\cdots,n_1)$') 
plot(real(lam_lambert),imag(lam_lambert),'o','DisplayName','Lambert W')
plot(real(lam_lf),imag(lam_lf),'.','DisplayName','LF')
xlim([-1 1]*Radius*2), ylim([-1 1]*Radius*2)
legend('show')
axis square
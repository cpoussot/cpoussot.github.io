function ZA = impedence_abs(alpha_fun,beta_fun,omega, porosity, pore_mean_size, pore_standard_dev)

% add 0 and infinite contributions
alpha           = alpha_fun([omega, porosity, pore_mean_size, pore_standard_dev]);
nu              = 1.5e-5;
tortuosity      = exp(4.*(pore_standard_dev.*log(2)).^2);
k0              = porosity.*pore_mean_size.^2./(8.*tortuosity) .* exp(-6*(pore_standard_dev.*log(2)).^2);    
M               = nu.*porosity ./ (k0.*tortuosity);
c_infty         = tortuosity./porosity;
alpha           = alpha(:) + tortuosity.*M./(porosity.*omega) + c_infty;
%
beta            = beta_fun([omega, porosity, pore_mean_size, pore_standard_dev]);
c_infty         = porosity;
beta            = beta + c_infty;
%
h   = .05;
Z0  = 348.84; %rho*c0 
c1  = 142800; %rho*gamma*P0  
c2  = 7.2857142857142855*1e-6; %rho/(gamma*P0) 
%Z   = -1i*sqrt(c1*alpha./beta).*cot(omega.*h.*sqrt(c2*alpha.*beta));
Z   = -1i*sqrt(c1*alpha./beta).*1./tan(imag(omega).*h.*sqrt(c2*alpha.*beta));
A   = 1-abs((Z-Z0)./(Z+Z0)).^2;
ZA  = [Z A];
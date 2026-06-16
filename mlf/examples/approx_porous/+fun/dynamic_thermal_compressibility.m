function beta_dyn = dynamic_thermal_compressibility(omega, porosity, pore_mean_size, pore_standard_dev)

    nu              = 1.5e-5;
    gamma           = 1.4;
    Pr              = 0.7;
    tortuosity      = exp(4*(pore_standard_dev.*log(2)).^2);
    k0_p            = porosity.*pore_mean_size.^2./(8.*tortuosity) .* exp(6.*(pore_standard_dev.*log(2)).^2);    
    thermal_length  = pore_mean_size .* exp(3/2.*(pore_standard_dev.*log(2)).^2);
    c_infty         = porosity;

    % JCAL
    M_p = nu*porosity ./ (k0_p.*Pr);
    N_p = nu*porosity ./ (k0_p.*Pr);
    L_p = nu .* porosity.^2 .* thermal_length.^2 ./ (4*k0_p.^2.*Pr);

    % % JCAPL
    % M_p = nu*porosity/(k0_p*Pr);
    % N_p = 2*nu/(thermal_length^2 *(static_thermal_tortuosity-1)*Pr);
    % L_p = nu/(thermal_length^2*(static_thermal_tortuosity-1)^2*Pr);

    % --- Dynamic viscous tortuosity
    beta_dyn = porosity.*gamma - porosity.*(gamma-1)./(1 ...
        + M_p./(omega) ...
        + N_p.*(sqrt(1+omega./L_p)-1)./(omega));

    beta_dyn = beta_dyn - c_infty;

end
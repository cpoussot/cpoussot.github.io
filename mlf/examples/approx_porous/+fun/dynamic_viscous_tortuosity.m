function alpha_dyn = dynamic_viscous_tortuosity(omega, porosity, pore_mean_size, pore_standard_dev)
    nu=1.5e-5;

    tortuosity = exp(4.*(pore_standard_dev.*log(2)).^2);
    k0 = porosity.*pore_mean_size.^2./(8.*tortuosity) .* exp(-6*(pore_standard_dev.*log(2)).^2);    
    viscous_length = pore_mean_size .* exp(-5/2.*(pore_standard_dev.*log(2)).^2);

    c_infty = tortuosity./porosity;

    % JCAL
    M = nu.*porosity ./ (k0.*tortuosity);
    N = nu.*porosity ./ (k0.*tortuosity);
    L = nu .* porosity.^2 .* viscous_length.^2 ./ (4.*k0.^2.*tortuosity.^2);

    % % JCAPL
    % M = nu*porosity/(k0*tortuosity);
    % N = 2*nu/(viscous_length^2 *(static_viscous_tortuosity/tortuosity-1));
    % L = nu/(viscous_length^2*(static_viscous_tortuosity/tortuosity-1)^2);
       
    % --- Dynamic viscous tortuosity
    alpha_dyn = tortuosity./porosity .*(1 ...
        + M./(omega) ...
        + N.*(sqrt(1+omega./L)-1)./(omega));

    alpha_dyn = alpha_dyn - tortuosity.*M./(porosity.*omega) - c_infty;
end

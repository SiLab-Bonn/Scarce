function Vbi = getDiffusionPotential(Neff, T)
    % Calculate the diffusion potential Vbi [V] as a function of the effective doping
    % concentration Neff [10^12 / cm^3] and te temperature [K]. 
    % Check/citetation of formulars needed!

    %constants
    Q_e = 1.602e-19;   %electron charge [C]
    K_b = 1.38e-23;   %Boltzmann Constant [J/K]
    
    N_i = 9.38.*10^7.*(T/300)^2.*exp(-6884/T);     % [10^12/cm^3] empirical fit at K = 300 range
    Vbi = K_b.*T./Q_e.*log(Neff.^2./N_i.^2);       % diffusion potential in thermal equilibrium [V]
 end %getTrappingTime  



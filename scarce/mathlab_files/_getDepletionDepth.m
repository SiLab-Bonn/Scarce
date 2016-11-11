function Xd = getDepletionDepth(Vbias, Neff, T)
    % Calculate depletion depth of silcon with an effective doping concentration of Neff [10^12 /cm^3]
    % at a temperature T [K] and with an reverse bias of Vbias [V]. 
    % Check/citetation of formulars needed!

    %constants
    Q_e = 1.602e-19;   %electron charge [C]
    K_b = 1.38e-23;   %Boltzmann Constant [J/K]
    Epsilon_0 = 8.85418781762e-12; %permittivity of vacuum [F/m]
    Epsilon_s = 1.04e-10; %permitivity of silicon [F/m]
    Epsilon_r = Epsilon_s/Epsilon_0; %relative permitivity of silicon [F/m]   
    
    Vbi= getDiffusionPotential(Neff, T);
    Xd = sqrt(2.*Epsilon_0.*Epsilon_r.*(Vbi+Vbias)./(Q_e.*Neff.*10.^6));    % depletion depth [m]
 end %getDepletionDepth


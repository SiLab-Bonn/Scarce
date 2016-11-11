function E_x = getElectricFieldPlanar(x, V_bias, N_eff, D)
    % Calculates the field E_x [V/um] in a planar sensor as a function of the position
    % x between the electrodes [um], the bias Voltage Vbias [V], the effective
    % doping concentration [Neff] and the sensor Width D [um]. 
    % The analytical function from the detector book p. 93 is used.

    V_dep = getDepletionVoltage(N_eff, D); %depletion voltage
    
    a = (V_bias - V_dep)./D;
    b = -2.*V_dep./(D.^2);
    E_x = (a-b.*x);
 end %getElectricFieldPlanar



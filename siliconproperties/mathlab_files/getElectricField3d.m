function [E_x E_y] = getElectricField3d(x, y, V_bias, N_eff, D, R)
    % Calculates the field E [V/um] in a 3d sensor as a function of the position
    % x between the electrodes [um], the bias Voltage Vbias [V], the effective
    % doping concentration [Neff], the electrode distance D [um] and radius R [um]. 
    % So far the same field like the weighting field is used --> space charge is ignored.

    V_dep = getDepletionVoltage(N_eff, D); %depletion voltage
    
    a = (V_bias - V_dep)./D;
    b = -2.*V_dep./(D.^2);
    [E_x E_y] = getWeightingField3d(x,y,D,R);
    E_x = E_x .* V_bias;
    E_y = E_y .* V_bias;
 end %getElectricField3d



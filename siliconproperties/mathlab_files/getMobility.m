function mu = getMobility(E, T, isElectron)
    % Calculate the mobility [cm^2/Vs] of charge carriers in silicon from 
    % the electrical fiel (E [V/cm]) the temperature (T [K]) and the charge
    % carrier type (isElectron [0/1] otherwise hole).
    % Formular derived from measured data of purity silicon and the 
    % corresponding fit function parameters, from:
    % C. Jacononi et al., Solid state electronics, 1977, vol 20., p. 87
    % 'A review of some charge transport properties of silicon'
    % the doping concentration is irrelevant for < 10^16/cm^3

    % constants
    if (isElectron==1) %electrons
        v_m = 1.53e9.*T.^(-0.87); %[cm/s]
        E_c = 1.01.*T.^(1.55); %[V/cm]
        beta = 2.57e-2.*T.^(0.66);  %[]
    else %holes
        v_m = 1.62e8.*T.^(-0.52); %[cm/s]
        E_c = 1.24.*T.^(1.68); %[V/cm]
        beta = 0.46.*T.^(0.17);  %[]
    end;
   
    mu = v_m./(E_c.* (1+(abs(E)./E_c).^beta).^(1/beta) );   
 end %getMobility  



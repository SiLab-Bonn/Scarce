function tr = getMeanFreePathWen(fluence)
    % Calculates the mean free path of electron holes in silicon.
    % Takes the data from A. Affolder, IEEE Nucl. Sc. Vol. 56, 3 June 2009
    % and calculates the mean free path from the CCE under the assumption
    % that the device is a pad detector (although it was a strip). 
    % Errors expected with this tough simplification...
    
    T = 300;    %temerature [K]
    N0 = 1.45; %doping concentration [10^12 /cm^3]
    isNtype = 1; %type of the sensor
    isOxygenated = 1; %oxygenated sensor
    D = 300; % detector width [um]
    S = 200; %pixel width[um]
    e_h = 1e5;  %deposited charge in electron hole pairs
    x_0 = 0;  %offset from the middle of the sensor pixel (x-direction) [um]
    y_0 = 150;  %position of the charge in the sensor depth/distance from readout electrode [um]
    V_dep = getDepletionVoltage(N0, D) % depletion voltage at starting doping concentration
    Vbias = -1.5*V_dep; %bias voltage [V];
    E = 1e6; % electric field [V/cm]

    tr_e = getTrappingTime(fluence, 1);
    tr_h = getTrappingTime(fluence, 0);

    tr = getTrappingTime(fluence, 1) .* getMobility(E, T, 1).*E;
 end %getMeanFreePathWen  



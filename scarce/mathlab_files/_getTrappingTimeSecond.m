function [tr tr_error] = getTrappingTime(fluence, isElectron)
    % Calculate the trapping time tr (e^-(tr) in ns) of charge carriers in 
    % silicon at a temperature of T = 263 K as a function of the 
    % fluence (in 10^12 Neq/cm^2). This was measured with a laser and planar 
    % silicon sensors with a fluence up to 2*10^14 Neq/cm^2. There was a 
    % linear behaviour between fluence and the effective trapping 
    % propability measured intepended of the silicon type (oxygenated or 
    % not and with different doping concentrations) from:
    % G. Kramberger et al., Nucl. Inst. And. Meth. A 476 (2002) 645-651
    % 'Determination of effective trapping times for electrons and holes 
    % in irradiated silicon'

    % constants
    if (isElectron==1) %electrons
%         beta = 4.2e-4;   %[cm^2/ns]
        beta = 2.98e-4;   %[cm^2/ns]
        beta_error = 0.3e-4;   %[cm^2/ns]
    else %holes
        beta = 3.57e-4;   %[cm^2/ns]
%         beta = 6.1e-4;   %[cm^2/ns]
        beta_error = 0.3e-4;   %[cm^2/ns]
    end;
   
    tr = 1./(fluence.*beta);
    tr_error = 1./(fluence.*beta^2).*beta_error;
 end %getTrappingTime  



function [tr tr_error] = getTrappingTime(fluence, isElectron)
    % Calculate the trapping time tr (e^-(tr) in ns) of charge carriers in 
    % silicon. From J. Weber Free Charge Carriers Trapping Properties in
    % Neutron-Irradiated DOFZ Silicon Pad Detectors

    % constants
    if (isElectron==1) %electrons
        beta = 4.2e-4;   %[cm^2/ns]
        beta_error = 0.3e-4;   %[cm^2/ns]
    else %holes
        beta = 6.1e-4;   %[cm^2/ns]
        beta_error = 0.3e-4;   %[cm^2/ns]
    end;
   
    tr = 1./(fluence.*beta);
    tr_error = 1./(fluence.*beta^2).*beta_error;
 end %getTrappingTime  



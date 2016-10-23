function [Neff Neff_error] = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated)
    % Calculate the effective acceptor concentration [10^12 cm^-3] of irradiated n and
    % p-type silicon with and without oxygen enriched as a function of the
    % fluence [10^12 cm^-2]. The data can be desribed by different line fits. 
    % The parameters were extracted from a plot taken from CERN/LHCC 
    % 2000-009 LEB Status Report/RD48 31 Dec. 1999 for the n-type silicon
    % and from RD50 Status Report 2006 CERN-LHCC-2007-005, p. 4-6 for p-type
    % silicon. Due to the difference in the data for different technologies
    % a rather large error on the propotionality factor of 10% is assumed.

    % constants
    if (isNtype==1) %p-type silicon
        beta_n = -5.5e-2;   %[cm^-1]       %slope before type inversion
        beta_n_fz = 2.83e-2;   %[cm^-1]     %slope after type inversion, standart silicon FZ
        beta_n_ofz = 0.94e-2;   %[cm^-1]    %slope after type inversion, oxygenated silicon FZ
        Phi_inv = -N0/beta_n; %inversion point N_eq [10^14 cm^-2]
        
        %after type inversion
        if(isOxygenated==1) %oxygenated silicon
            Neff = beta_n_ofz.*(fluence-Phi_inv);
        else %standart silicon
            Neff = beta_n_fz.*(fluence-Phi_inv);
        end;      
        
        %before type inversion
        Neff(fluence<Phi_inv) = N0+beta_n.*fluence(fluence<Phi_inv);

    else %n-type silicons
        beta_n_fz = 2.3e-2;   %[cm^-1]
        beta_n_ofz = 2.7e-2;   %[cm^-1]
        if(isOxygenated==1) %oxygenated silicon
            Neff = N0+ beta_n_ofz.*fluence;
        else %standart silicon
            Neff = N0+ beta_n_fz.*fluence;
        end;
    end;
    Neff_error = Neff.*0.1;
 end %getEffAcceptorConentration



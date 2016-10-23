function V_dep = getDepletionVoltage(Neff, Electroded)
    %this function returns the needed depletion Voltage [V] as a function
    %of the effective doping concentration Neff [10^12 /cm^-3] and the
    %distance between electrodes [um]. Formular is standard and for example
    %used in G. Kramberger et al., Nucl. Inst. And. Meth. A 476 (2002) 645-651
    % 'Determination of effective trapping times for electrons and holes 
    % in irradiated silicon'.

    %constants
    Epsilon_0 = 8.85418781762e-12; %permittivity of vacuum [F/m]
    Epsilon_r = 11.75; %relative permittivity of silicon [F/m]
    Q_e = 1.602e-19;   %electron charge [C]
    
    V_dep = Q_e.*Neff./(Epsilon_0.*Epsilon_r).*Electroded.^2./2*1e6;
end %getDepletionVoltage
clear;
%constants
Detector_density = 2.3290;  %density of detector material silicon [g cm^-3]
Doping_concentration = 1e12;    %effective doping concentration of silicon bulk(Nd-Na)
N_a = 1*10^15;    %acceptor doping concentration [1/cm^3]
N_d = 1*10^15;    %donator doping concentration [1/cm^3]
T = 300;  %temperature [K]
K_b = 1.38e-23;   %Boltzmann Constant [J/K]
Q_e = 1.602e-19;   %electron charge [C]
Epsilon_0 = 8.85418781762e-12; %permittivity of vacuum [F/m]
Epsilon_s = 1.04e-10; %permitivity of silicon [F/m]
Epsilon_r = Epsilon_s/Epsilon_0; %relative permitivity of silicon [F/m]
%the following need to be checked
E_g = 1.17 - 4.73*10^(-4)*T^2./(T+636);     %band gab energy of silicon [eV]
N_c = 4.82*10^(15)*6*(0.36)^(3/2)*T^(3/2);  % effective density of states in the conduction band [cm^-3]
N_v = 3.5*10^(15)*T^(3/2);                  % effective density of states in the valence band [cm^-3]
N_i =(N_c.*N_v )^(1/2).*exp(-E_g.*Q_e./(2*K_b*T));  % instrinsic carrier concentration of silicon
N_i = 9.38.*10^(19).*(T/300)^2.*exp(-6884/T);       % empirical fit at K = 300 range
V_bi = K_b.*T./Q_e.*log(N_a.*N_d./N_i^2);           % diffusion potential in thermal equilibrium [V]
X_d = sqrt(2.*Epsilon_0.*Epsilon_s.*V_bi./(Q_e.*N_d./10.^6))  % depletion depth without bias [m]
bias = 0:0.1:100;    % bias in volt
X_d_bias = sqrt(2.*Epsilon_0.*Epsilon_r.*(V_bi+bias)./(Q_e.*N_d./10.^6));    % depletion depth [m]

plot(bias,sqrt(2.*Epsilon_0.*Epsilon_r.*(V_bi+bias)./(Q_e.*1e11/10.^6)), 'LineWidth', 2, 'color', 'red');
hold on;
plot(bias,sqrt(2.*Epsilon_0.*Epsilon_r.*(V_bi+bias)./(Q_e.*3e11/10.^6)), 'LineWidth', 2, 'color', 'blue', 'LineStyle', '--');
plot(bias,sqrt(2.*Epsilon_0.*Epsilon_r.*(V_bi+bias)./(Q_e.*5e11/10.^6)), 'LineWidth', 2, 'color', 'green', 'LineStyle', '-.');
plot(bias,sqrt(2.*Epsilon_0.*Epsilon_r.*(V_bi+bias)./(Q_e.*1e12/10.^6)), 'LineWidth', 2, 'color', 'black', 'LineStyle', ':');
hold off;
title('Depletion depth for different doping concentrations ', 'FontWeight','bold','FontSize', 10);
xlabel('depletion voltage [V]', 'FontWeight','bold');
ylabel('depletion depth [um]', 'FontWeight','bold');
legend('1*10^{11}/cm^{3}','3*10^{11}/cm^{3}','5*10^{11}/cm^{3}','1*10^{12}/cm^{3}');
set(gcf, 'Color', [1 1 1]);
%export_fig('Depletion_Depth.pdf');
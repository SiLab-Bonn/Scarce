% be carefull constant max. saturation velocity assumed!

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

fluence = linspace(0,25e3,1e3);
tr_e = getTrappingTime(fluence, 1);
tr_h = getTrappingTime(fluence, 0);

lambda_e = getTrappingTime(fluence, 1) .* getMobility(E, T, 1).*E;
getMobility(E, T, 1).*E
getTrappingTime(5000, 1)

%semilogy(fluence,tr_e, 'LineWidth', 1.2, 'color', 'blue');
plot(fluence./10^3,lambda_e.*1e-4, 'LineWidth', 1.2, 'color', 'blue');
% hold on;
% %semilogy(fluence,tr_h, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');
% loglog(fluence,tr_h, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');
% hold off;
title('Charge carrier mean free path in irradiated silicon with different fluence', 'FontWeight','bold','FontSize', 10);
xlabel('fluence [10^{15} N_{eq}/cm^2]', 'FontWeight','bold');
ylabel('mean free path [um]', 'FontWeight','bold');
set(gcf, 'Color', [1 1 1]);
ylim([0 1e3]);
hleg1 = legend('electrons','holes', 'Location', 'northeast');
grid on;
set(gca, 'GridLineStyle', '-');
% export_fig('charge_carrier_mean_free_path.pdf');




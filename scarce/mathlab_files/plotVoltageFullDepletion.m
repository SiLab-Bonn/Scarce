%constants
clear;
N0 = 1;  %doping concentration [10^12 1/cm^3]
isNtype = 1;    %n-Type sensor
isOxygenated = 1;   %oxygenated sensor
d = 200;    %detector thickness [um]

% fluence = [0:1:350];
fluence = logspace(0,4,1000);
Neff = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated);
Vfd = getDepletionVoltage(Neff, d);

loglog(fluence,Vfd, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');

hold on;
% loglog(fluence,Neff_p_OFZ, 'LineWidth', 1.2, 'color', 'black');
% loglog(fluence,Neff_n_FZ, 'LineWidth', 1.2, 'color', 'red', 'LineStyle', '-.');
% loglog(fluence,Neff_n_OFZ, 'LineWidth', 1.2, 'color', 'green', 'LineStyle', ':');
% plot(fluence,Neff_p_OFZ, 'LineWidth', 1.2, 'color', 'black');
% plot(fluence,Neff_n_FZ, 'LineWidth', 1.2, 'color', 'red', 'LineStyle', '-.');
% plot(fluence,Neff_n_OFZ, 'LineWidth', 1.2, 'color', 'green', 'LineStyle', ':');
hold off;
title('Depletion voltage for oxygenated n-type 1*10^{12}/cm^3 planar silicon sensor', 'FontWeight','bold','FontSize', 10);
xlabel('fluence [10^{12} N_{eq}/cm^2]', 'FontWeight','bold');
ylabel('V_{D} [V]', 'FontWeight','bold');
set(gcf, 'Color', [1 1 1]);
ylim([1e-1 1e4]);
%hleg1 = legend('p-type FZ silicon','p-type oxigenated FZ silicon', 'n-type FZ silicon','n-type oxigenated FZ silicon', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', '-');
% export_fig('depletion_voltage.pdf');




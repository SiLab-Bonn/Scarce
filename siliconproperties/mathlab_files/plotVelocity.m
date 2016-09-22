E = logspace(4,6,1000);
mu0_e = getMobility(E, 245, 0);
mu1_e = getMobility(E, 300, 0);
mu2_e = getMobility(E, 370, 0);
mu3_e = getMobility(E, 430, 0);
mu0_h = getMobility(E, 245, 1);
mu1_h = getMobility(E, 300, 1);
mu2_h = getMobility(E, 370, 1);
mu3_h = getMobility(E, 430, 1);
semilogx(E,mu0_e.*E, 'LineWidth', 1.2, 'color', 'black');
hold on;
semilogx(E,mu1_e.*E, 'LineWidth', 1.2, 'color', 'red');
semilogx(E,mu2_e.*E, 'LineWidth', 1.2, 'color', 'blue');
semilogx(E,mu3_e.*E, 'LineWidth', 1.2, 'color', 'green');
semilogx(E,mu0_h.*E, 'LineWidth', 1.2, 'color', 'black', 'LineStyle', '--');
semilogx(E,mu1_h.*E, 'LineWidth', 1.2, 'color', 'red', 'LineStyle', '--');
semilogx(E,mu2_h.*E, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');
semilogx(E,mu3_h.*E, 'LineWidth', 1.2, 'color', 'green', 'LineStyle', '--');
hold off;
title('Charge carrier velocity in high purity silicon ', 'FontWeight','bold','FontSize', 10);
xlabel('electric field [V/cm]', 'FontWeight','bold');
ylabel('electron/-hole velocity [cm/s]', 'FontWeight','bold');
set(gcf, 'Color', [1 1 1]);
% ylim([1e4 2e7]);
ylim([1e4 1.4e7]);
hleg1 = legend('T = 245 K','T = 300 K','T = 370 K','T = 430 K', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', '-');
% export_fig('charge_carrier_velocity.pdf');
export_fig('charge_carrier_velocity_zoom.pdf');
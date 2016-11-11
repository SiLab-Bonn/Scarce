fluence = 0:10:10000;
% fluence = logspace(2,16,1000);
tr_e = getTrappingTime(fluence, 1);
tr_h = getTrappingTime(fluence, 0);
semilogy(fluence,tr_e, 'LineWidth', 1.2, 'color', 'blue');
% loglog(fluence,tr_e, 'LineWidth', 1.2, 'color', 'blue');
hold on;
semilogy(fluence,tr_h, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');
% loglog(fluence,tr_h, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');
hold off;
title('Charge carrier trapping time in irradiated silicon with different fluence', 'FontWeight','bold','FontSize', 10);
xlabel('fluence [10^{12} N_{eq}/cm^2]', 'FontWeight','bold');
ylabel('trapping time [ns]', 'FontWeight','bold');
set(gcf, 'Color', [1 1 1]);
% ylim([1e-1 1e4]);
hleg1 = legend('electrons','holes', 'Location', 'northeast');
grid on;
set(gca, 'GridLineStyle', '-');
% export_fig('charge_carrier_trapping_time.pdf');




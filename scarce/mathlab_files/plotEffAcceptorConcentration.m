% fluence = [0:1:350];
clear;
fluence = logspace(0,4,1000);
Neff_p_FZ = getEffAcceptorConentration(fluence, 1.8, 0, 0);
Neff_p_OFZ = getEffAcceptorConentration(fluence, 1.8, 0, 1);
Neff_n_FZ = getEffAcceptorConentration(fluence, 1.8, 1, 0);
Neff_n_OFZ = getEffAcceptorConentration(fluence, 1.8, 1, 1);
loglog(fluence,Neff_p_FZ, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');
% plot(fluence,Neff_p_FZ, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');
hold on;
loglog(fluence,Neff_p_OFZ, 'LineWidth', 1.2, 'color', 'black');
loglog(fluence,Neff_n_FZ, 'LineWidth', 1.2, 'color', 'red', 'LineStyle', '-.');
loglog(fluence,Neff_n_OFZ, 'LineWidth', 1.2, 'color', 'green', 'LineStyle', ':');
% plot(fluence,Neff_p_OFZ, 'LineWidth', 1.2, 'color', 'black');
% plot(fluence,Neff_n_FZ, 'LineWidth', 1.2, 'color', 'red', 'LineStyle', '-.');
% plot(fluence,Neff_n_OFZ, 'LineWidth', 1.2, 'color', 'green', 'LineStyle', ':');
hold off;
title('Effective acceptor concentration for silicon after irradiation', 'FontWeight','bold','FontSize', 10);
xlabel('fluence [10^{12} N_{eq}/cm^2]', 'FontWeight','bold');
ylabel('N_{eff} [10^{12} cm^{-3}]', 'FontWeight','bold');
set(gcf, 'Color', [1 1 1]);
ylim([1e-2 1e3]);
hleg1 = legend('p-type FZ silicon','p-type oxigenated FZ silicon', 'n-type FZ silicon','n-type oxigenated FZ silicon', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', '-');
%export_fig('effective_acceptor_concentration_zoom.pdf');
% export_fig('effective_acceptor_concentration.pdf');




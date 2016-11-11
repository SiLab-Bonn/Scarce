%constants
Neff = 0.1;  %doping concentration [10^12 1/cm^3]
Vbias = 0:0.1:100; %depletion voltage
T = 300;    %temperature [K]

plot(Vbias,getDepletionDepth(Vbias, Neff, T), 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');

hold on;
plot(Vbias,getDepletionDepth(Vbias, 5*Neff, T), 'LineWidth', 1.2, 'color', 'black');
plot(Vbias,getDepletionDepth(Vbias, 10*Neff, T), 'LineWidth', 1.2, 'color', 'red', 'LineStyle', '-.');
plot(Vbias,getDepletionDepth(Vbias, 15*Neff, T), 'LineWidth', 1.2, 'color', 'green', 'LineStyle', ':');
hold off;
title('Depletion depth for different doping concentrations ', 'FontWeight','bold','FontSize', 10);
xlabel('depletion voltage [V]', 'FontWeight','bold');
ylabel('depletion depth [um]', 'FontWeight','bold');
set(gcf, 'Color', [1 1 1]);
legend('1*10^{11}/cm^{3}','5*10^{11}/cm^{3}','1*10^{12}/cm^{3}','5*10^{12}/cm^{3}', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', '-','XMinorTick','on','YMinorTick','on');
export_fig('Depletion_Depth.pdf');




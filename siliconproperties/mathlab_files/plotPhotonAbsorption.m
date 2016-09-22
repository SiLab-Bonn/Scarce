clear;
%constants
detector_density = 2.3290;  %density of detector material [g cm^-3]

fileID = fopen('AttenuationData.txt');
%file format:
%Photon   	Coherent	Photoel.	Nuclear 	
%Energy   	Scatter.	Absorb. 	Pr. Prd.
data = textscan(fileID, '%f %f %f %f %f');
loglog(data{1}, data{2}.*detector_density, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '--');
hold on;
loglog(data{1}, data{3}.*detector_density, 'LineWidth', 1.2, 'color', 'green', 'LineStyle', ':');
loglog(data{1}, data{4}.*detector_density, 'LineWidth', 1.2, 'color', 'red', 'LineStyle', '-.');
loglog(data{1}, data{5}.*detector_density, 'LineWidth', 1.2, 'color', 'black');
hold off;
title('Attenuation for photons in silicon (Z=14)', 'FontWeight','bold');
xlabel('photon energy [MeV]', 'FontWeight','bold');
ylabel('attenuation [1/cm]', 'FontWeight','bold');
xlim([0.001 10000]);
ylim([0.001 1e4]);
set(gcf, 'Color', [1 1 1]);
hleg1 = legend('compton effect','photo effect','pair production','total');
grid on
set(gca, 'GridLineStyle', '-');
export_fig('photon_absoption.pdf');
%close(fileID);
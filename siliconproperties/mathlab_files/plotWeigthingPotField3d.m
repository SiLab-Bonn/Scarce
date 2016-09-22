% constants
S = 70; %electrode distance [um]
R = 2.5; %electrode radius [um]
x = [-S:0.001:S];  %position in the sensor depth [um]
y = zeros(1, size(x,2));    %get field/potential on the shortest line between the electrodes (y = 0)

[E_x_w E_y_w] = getWeightingField3d(x,y,S,R);
[AX1,H1,H2] = plotyy(x,getWeightingPotential3d(x,y,S,R),x,E_x_w);
hold on;
[E_x_w E_y_w] = getWeightingField3d(x,y,S,R*4);
[AX2,HH1,HH2] = plotyy(x,getWeightingPotential3d(x,y,S,R*4),x, E_x_w);
% [E_x_w E_y_w] = getWeightingFieldPlanar(x,y,D,S*16);
% [AX3,HHH1,HHH2] = plotyy(y,getWeightingPotentialPlanar(0,y,D,S*16),y,-E_y_w);
hold off;
title_str = sprintf('Weighting potential and field, one 3d pixel, distance %d um', S);
title(title_str, 'FontWeight','bold','FontSize', 10);
xlabel('position [um]', 'FontWeight','bold');
set(get(AX1(1), 'ylabel'), 'String','weighting potential [V]'); %, 'FontWeight','bold');
set(get(AX1(2), 'ylabel'), 'String','weighting field [V/um]');
set(H1, 'LineWidth', 1.2);
set(H2, 'LineWidth', 1.2);
set(HH1, 'LineWidth', 1.2, 'LineStyle', '--');
set(HH2, 'LineWidth', 1.2, 'LineStyle', '--');
% set(HHH1, 'LineWidth', 1.2, 'LineStyle', '-.');
% set(HHH2, 'LineWidth', 1.2, 'LineStyle', '-.');
set(AX1(2), 'ylim', [-0.10 0.10]);
set(AX2(2), 'ylim', [-0.10 0.10]);
set(AX1, 'xlim', [-100 100]);
set(AX2, 'xlim', [-100 100]);
% set(AX3(2), 'ylim', [0 0.06]);
set(AX1(1),'YTick',[0:0.1:1]);
set(AX1(2),'YTick',[-0.10:0.01:0.10]);
set(AX2(2),'YTick',[-0.10:0.01:0.10]);
% set(AX3(2),'YTick',[0:0.004:0.06]);
set(gcf, 'Color', [1 1 1]);
%ylim([1e-2 1e3]);
hleg1 = legend('\Phi_{w}, |E|_{w}, col. radius 2.5 um','\Phi_{w}, |E|_{w}, col. radius 10 um', 'Location', 'southwest');
grid on;
set(gca, 'GridLineStyle', '-');
export_fig('WeightingField_3d.pdf');




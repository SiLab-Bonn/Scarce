% constant
D = 300; %detector width [um]
S = 25; %electrode size [um]
y = [0:0.1:D];  %position in the sensor depth [um]

x = zeros(1, size(y,2));    %get field/potential in the middle of electrode (x = 0)

[E_x_w E_y_w] = getWeightingFieldPlanar(x,y,D,S);
[AX1,H1,H2] = plotyy(y,getWeightingPotentialPlanar(0,y,D,S),y,-E_y_w);
hold on;
[E_x_w E_y_w] = getWeightingFieldPlanar(x,y,D,S*4);
[AX2,HH1,HH2] = plotyy(y,getWeightingPotentialPlanar(0,y,D,S*4),y, -E_y_w);
[E_x_w E_y_w] = getWeightingFieldPlanar(x,y,D,S*16);
[AX3,HHH1,HHH2] = plotyy(y,getWeightingPotentialPlanar(0,y,D,S*16),y,-E_y_w);
hold off;
title_str = sprintf('Weighting potential and field, one planar pixel, thickness %d um, pixel size %d um', D, S);
title(title_str, 'FontWeight','bold','FontSize', 10);
xlabel('position [um]', 'FontWeight','bold');
set(get(AX1(1), 'ylabel'), 'String','potential [V]'); %, 'FontWeight','bold');
set(get(AX1(2), 'ylabel'), 'String','electrical field [V/um]');
set(H1, 'LineWidth', 1.2);
set(H2, 'LineWidth', 1.2);
set(HH1, 'LineWidth', 1.2, 'LineStyle', '--');
set(HH2, 'LineWidth', 1.2, 'LineStyle', '--');
set(HHH1, 'LineWidth', 1.2, 'LineStyle', '-.');
set(HHH2, 'LineWidth', 1.2, 'LineStyle', '-.');
set(AX1(2), 'ylim', [0 0.06]);
set(AX2(2), 'ylim', [0 0.06]);
set(AX3(2), 'ylim', [0 0.06]);
set(AX1(1),'YTick',[0:0.1:1]);
set(AX1(2),'YTick',[0:0.004:0.06]);
set(AX2(2),'YTick',[0:0.004:0.06]);
set(AX3(2),'YTick',[0:0.004:0.06]);
set(gcf, 'Color', [1 1 1]);
%ylim([1e-2 1e3]);
hleg1 = legend('\Phi_{w}, |E|_{w}, pixel width 25 um','\Phi_{w}, |E|_{w}, pixel width 100 um', '\Phi_{w}, |E|_{w}, pixel width 400 um', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', '-');
export_fig('WeightingField_planar.pdf');




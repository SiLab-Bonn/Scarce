%constants
N_eff = 1;  %doping concentration [10^12 1/cm^3]
D = 200; %   sensor area in y plot [um]
S = 50;    %sensor columns distance [um]
R = 5; %sensor columns radius [um]
V_bias = -getDepletionVoltage(N_eff, S); %bias voltage
x = -100:1:100; %position in sensor [um]
y = 0:1:200; %position in sensor [um]
[xx,yy] = meshgrid(x,y);

[E_x E_y] = getElectricField3d(xx, yy-D/2, V_bias, N_eff, S, R);
E_x(E_x>1) = 1;
E_y(E_y>1) = 1;
E_x(E_x<-1) = NaN;
E_y(E_y<-1) = NaN;


% E_y_2 = getElectricFieldPlanar(y, 0.5.*V_bias, N_eff, D);
% E_y_2(E_y_2>0) = 0;
% 
% E_y_3 = getElectricFieldPlanar(y, 2.*V_bias, N_eff, D);
% E_y_3(E_y_3>0) = 0;
% 
% E_y_4 = getElectricFieldPlanar(y, V_bias, 0.*N_eff, D);
% E_y_4(E_y_4>0) = 0;
% 
% E_y_5 = getElectricFieldPlanar(y, V_bias, 0.5.*N_eff, D);
% E_y_5(E_y_5>0) = 0;
% 
% E_y_6 = getElectricFieldPlanar(y, 10*V_bias, 10.5*N_eff, D);
% E_y_6(E_y_6>0) = 0;

plot(x, E_x(floor(size(E_x,1)/2)+1,:), 'LineWidth', 1.2, 'color', 'black', 'LineStyle', '-');

% hold on;
% plot(y, E_y_2, 'LineWidth', 1.2, 'color', 'black', 'LineStyle', ':');
% plot(y, E_y_3, 'LineWidth', 1.2, 'color', 'black', 'LineStyle', '-.');
% plot(y, E_y_4, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '-');
% plot(y, E_y_5, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', ':');
% plot(y, E_y_6, 'LineWidth', 1.2, 'color', 'blue', 'LineStyle', '-.');
% hold off;
title('Electric field between two columns in a 3d silicon sensor ', 'FontWeight','bold','FontSize', 10);
xlabel('position [um]', 'FontWeight','bold');
ylabel('electric field along x [V/um]', 'FontWeight','bold');
set(gcf, 'Color', [1 1 1]);
% legend('1*10^{12}/cm^{3}, V_{bias}=V_{dep}','1*10^{12}/cm^{3}, V_{bias}=0.5 V_{dep}','1*10^{12}/cm^{3}, V_{bias}=2 V_{dep}','0*10^{12}/cm^{3}, V_{bias}=V_{dep}' ,'0.5*10^{12}/cm^{3}, V_{bias}=V_{dep}','1.5*10^{12}/cm^{3}, V_{bias}=V_{dep}','Location', 'southeast');
grid on;
set(gca, 'GridLineStyle', '-','XMinorTick','on','YMinorTick','on');
%export_fig('Electric_Field_Planar.pdf');




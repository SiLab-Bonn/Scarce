clear;
plotCC3d;    
hold on;
A = importdata('CCE_3d.csv');
plot(A(:,1),A(:,3),'.','color','red','MarkerSize',10);
hold on;
plotCCplanar;    
hold on;
A = importdata('CCE_planar.csv');
plot(A(:,1),A(:,3),'.','color','blue','MarkerSize',10);
hold on;

fluence = 0:1e3:2e4;    %radiation dose [10^12 Neq/cm^3]

T = 248;    %temperature [K]
N0 = 1.00; %doping concentration [10^12 /cm^3]
isNtype = 1; %type of the sensor
isOxygenated = 1; %oxygenated sensor
D = 70; % detector width [um]
S = 100; %pixel width[um]
Qin = 4.37e3;  %deposited charge in electron hole pairs
xtrack = 0;  %offset from the middle of the sensor pixel (x-direction) [um]
Neff = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated);
Vbias_1 = -80; %bias voltage [V];
dt = 1e-3;   %time step of simulation [ns]
t = 0:dt:6; %time points of simulation [ns]
N = 50;

%plot CCE as a function of the fluence
CC_1=zeros(1,size(fluence,2));
CCE_2=zeros(1,size(fluence,2));
lambda_1=zeros(1,size(fluence,2));
lambda_2=zeros(1,size(fluence,2));
iter = 1;

for(ifluence = fluence)
	Neff = getEffAcceptorConentration(ifluence, N0, isNtype, isOxygenated);
	Vbias_1 = -1.1*getDepletionVoltage(Neff, D)-20;
	Q_tot_1 = getSignalPlanarSensorSimple(xtrack,Qin,D,S,N0,isNtype,isOxygenated,Vbias_1,ifluence,T, t,dt, N);
% 	Vbias_2 = -getDepletionVoltage(Neff, D)-20;
% 	Q_tot_2 = getSignalPlanarSensor(xtrack,Qin,D,S,N0,isNtype,isOxygenated,Vbias_2,ifluence,T, t,dt, N);
	CC_1(iter)=abs(min(Q_tot_1))./1e3;
% 	CC_2(iter)=abs(min(Q_tot_2))./1e3;
	iter = iter + 1;
	ifluence
    Neff
end
        
plot(fluence./1e3, CC_1, 'COLOR', 'green','LineWidth', 1.2, 'LineStyle', '-');
% plot(fluence./1e3, CC_2, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '--');

hold off;
title('Collected charge for planar and 3d silicon pixel sensors, 300 thickness', 'FontWeight','bold','FontSize', 10);
xlabel('fluence [10^{15} N_{eq}/cm^2]', 'FontWeight','bold');
ylabel('collected charge [ke]', 'FontWeight','bold');
legend('3d calculated','3d measured','planar calculated','planar calc. enhanced','planar measured','planar calc. with 70 um thickness','Location', 'northeast');
% legend('3d calculated','3d measured','planar calculated','planar measured','Location', 'northeast');
% legend('3d measured','planar measured','Location', 'northeast');
grid on;
set(gca, 'GridLineStyle', '-');
ylim([0 24]);
set(gcf, 'Position', [100 100 100+600 100+350], 'Color', [1 1 1]);
% export_fig('CC_comparison_measured.pdf');
% export_fig('CC_comparison.pdf');
hold off;

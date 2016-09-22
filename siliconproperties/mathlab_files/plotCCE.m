clear;
plotCCE3d;    
hold on;
A = importdata('CCE_3d.csv');
plot(A(:,1),A(:,2),'.','color','red','MarkerSize',10);
hold on;
plotCCEplanar;    
hold on;
A = importdata('CCE_planar.csv');
plot(A(:,1),A(:,2),'.','color','blue','MarkerSize',10);

fluence = 0:1e3:2e4;    %radiation dose [10^12 Neq/cm^3]

T = 248;    %temperature [K]
N0 = 1.00; %doping concentration [10^12 /cm^3]
isNtype = 1; %type of the sensor
isOxygenated = 1; %oxygenated sensor
D = 70; % detector width [um]
S = 100; %pixel width[um]
Qin = 1e5;  %deposited charge in electron hole pairs
xtrack = 0;  %offset from the middle of the sensor pixel (x-direction) [um]
Neff = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated);
Vbias_1 = -80; %bias voltage [V];
dt = 1e-3;   %time step of simulation [ns]
t = 0:dt:6; %time points of simulation [ns]
N = 50;

%plot CCE as a function of the fluence
CCE_1=zeros(1,size(fluence,2));
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
	CCE_1(iter)=abs(min(Q_tot_1)/Qin)*100;
% 	CCE_2(iter)=abs(min(Q_tot_2)/Qin)*100;
	iter = iter + 1;
% 	ifluence
%     Neff
end
        
plot(fluence./1e3, CCE_1, 'COLOR', 'green','LineWidth', 1.2, 'LineStyle', '--');

hold off;
title('Charge collection efficiencies for planar/3d silicon pixel sensors, 300/70 um electrode distance', 'FontWeight','bold','FontSize', 10);
xlabel('fluence [10^{15} N_{eq}/cm^2]', 'FontWeight','bold');
ylabel('CCE [%]', 'FontWeight','bold');
legend('3d calculated','3d measured','planar measured','planar measured','planar calculated 70 um','Location', 'northeast');
% legend('3d sensor','planar sensor','Location', 'northeast');
grid on;
set(gca, 'GridLineStyle', '-');
ylim([0 100]);
set(gcf, 'Position', [100 100 100+600 100+350], 'Color', [1 1 1]);
% export_fig('CCE_comparison.pdf');
% export_fig('CCE_comparison_measured.pdf');
hold off;

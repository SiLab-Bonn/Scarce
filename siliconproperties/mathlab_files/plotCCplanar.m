%constants
fluence = 0:1e3:2e4;    %radiation dose [10^12 Neq/cm^3]

T = 248;    %temperature [K]
N0 = 1.00; %doping concentration [10^12 /cm^3]
isNtype = 1; %type of the sensor
isOxygenated = 1; %oxygenated sensor
D = 300; % detector width [um]
S = 100; %pixel width[um]
Qin = 23e3;  %deposited charge in electron hole pairs
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
	Vbias_2 = -getDepletionVoltage(Neff, D)-20;
	Q_tot_2 = getSignalPlanarSensor(xtrack,Qin,D,S,N0,isNtype,isOxygenated,Vbias_2,ifluence,T, t,dt, N);
	CC_1(iter)=abs(min(Q_tot_1))./1e3;
	CC_2(iter)=abs(min(Q_tot_2))./1e3;
	iter = iter + 1;
	ifluence
    Neff
end
        
plot(fluence./1e3, CC_1, 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '-');
hold on;
plot(fluence./1e3, CC_2, 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '--');
% hold off;
title('Collected charge of a planar 300 um pixel detector with 100 um pixel', 'FontWeight','bold','FontSize', 10);
xlabel('fluence [10^{15} N_{eq}/cm^2]', 'FontWeight','bold');
ylabel('collected charge [ke]', 'FontWeight','bold');
% legend('V_{bias} = 1.1*V_{dep}-20','Location', 'northeast');
grid on;
set(gcf, 'Color', [1 1 1]);
set(gca, 'GridLineStyle', '-');
ylim([0 24]);
hold off;

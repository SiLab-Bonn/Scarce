%constants
clear;
fluence = 0:1e3:2e4;    %radiation dose [10^12 Neq/cm^3]

T = 248;    %temperature [K]
N0 = 1.00; %doping concentration [10^12 /cm^3]
isNtype = 1; %type of the sensor
isOxygenated = 1; %oxygenated sensor
D = 300; % detector width [um]
S = 70; %electrode distance [um]
R = 5; %column radius [um]
Qin = 23e3;  %deposited charge in electron hole pairs
% xtrack = 0;  %offset from the middle of the sensor pixel (x-direction) [um]
% ytrack = D/2;  %offset from the middle of the sensor pixel (y-direction) [um]
xtrack = linspace(-S/2+R,S/2-R,10);  %offset from the middle of the sensor pixel (x-direction) [um]
xtrack = repmat(xtrack,[size(xtrack,2) 1]);
xtrack = xtrack(:)';
ytrack = linspace(D/2-15,D/2+15,floor(sqrt(size(xtrack,2))));  %offset from the bottom of the sensor pixel (y-direction) [um]
ytrack = repmat(ytrack,1,size(ytrack,2));
Neff = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated);
Vbias_1 = -1.5*getDepletionVoltage(Neff, D); %bias voltage [V];
dt = 1e-3;   %time step of simulation [ns]
t = 0:dt:1.0; %time points of simulation [ns]
%     t = 0;
Resolution = 200;

%plot CCE as a function of the fluence
CCE_1=zeros(1,size(fluence,2));
% CCE_2=zeros(1,size(fluence,2));
lambda_1=zeros(1,size(fluence,2));
% lambda_2=zeros(1,size(fluence,2));
iter = 1;

for(ifluence = fluence)
	Neff = getEffAcceptorConentration(ifluence, N0, isNtype, isOxygenated);
	Vbias_1 = 0;
	[Q_ind_tot Q_ind_e_vec Q_ind_h_vec] = getSignal3dSensor(xtrack,ytrack,Qin,D,S,R,N0,isNtype,isOxygenated,Vbias_1,ifluence,T, t,dt);
% 	Vbias_2 = -getDepletionVoltage(Neff, D)-20;
% 	Q_tot_2 = getSignal3dSensor(xtrack,ytrack,Qin,D,S,R,N0,isNtype,isOxygenated,Vbias_2,ifluence,T, t,dt);
	CCE_1(iter)=abs(min(mean(Q_ind_e_vec(:,:)+Q_ind_h_vec(:,:),1)))./1e3;
% 	CCE_2(iter)=abs(min(Q_tot_2)/Qin)*100;
	iter = iter + 1;
	ifluence
end
        
% plot(t, Q_tot_1./1e5, 'COLOR', 'blue','LineWidth', 2, 'LineStyle', '-');   

plot(fluence./1e3, CCE_1, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '-');
hold on;
% plot(fluence, CCE_2, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '-');
% hold off;
title_str = sprintf('Collected charge of a 3d 300 um pixel detector, %1.0f um electrode distance, %1.0f um radius', S, R);
title(title_str, 'FontWeight','bold','FontSize', 10);
xlabel('fluence [10^{15} N_{eq}/cm^2]', 'FontWeight','bold');
ylabel('collected charge [ke]', 'FontWeight','bold');
% legend('V_{bias} = 1.1*V_{dep}-20','Location', 'northeast');
grid on;
set(gcf, 'Color', [1 1 1]);
set(gca, 'GridLineStyle', '-');
ylim([0 24]);
hold off;

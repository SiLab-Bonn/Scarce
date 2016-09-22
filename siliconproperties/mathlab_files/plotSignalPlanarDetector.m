%constants
    clear;
    T = 248;    %temperature [K]
    fluence = 0e4;    %radiation dose [10^12 Neq/cm^3]
    N0 = 1.00; %doping concentration [10^12 /cm^3]
    isNtype = 1; %type of the sensor
    isOxygenated = 1; %oxygenated sensor
    D = 300; % detector width [um]
    S = 100; %pixel width[um]
    Qin = 3e5;  %deposited charge in electron hole pairs
    xtrack = 0;  %offset from the middle of the sensor pixel (x-direction) [um]
    Neff = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated);
    Vbias_1 = -1.5*getDepletionVoltage(Neff, D); %bias voltage [V];
    dt = 1e-3;   %time step of simulation [ns]
    t = 0:dt:5; %time points of simulation [ns]
    N = 3;
    
    %plot CCE as a function of the fluence
%         fluence = 0:1e3:1e4;    %radiation dose [10^12 Neq/cm^3]
%         CCE_1=zeros(1,size(fluence,2));
%         CCE_2=zeros(1,size(fluence,2));
%         lambda_1=zeros(1,size(fluence,2));
%         lambda_2=zeros(1,size(fluence,2));
%         iter = 1;
% 
%         for(ifluence = fluence)
%             Neff = getEffAcceptorConentration(ifluence, N0, isNtype, isOxygenated);
%             Vbias_1 = -1.1*getDepletionVoltage(Neff, D);
%             Q_tot_1 = getSignalPlanarSensor(xtrack,Qin,D,S,N0,isNtype,isOxygenated,Vbias_1,ifluence,T, t,dt, N);
%             Vbias_2 = -800;
%             Q_tot_2 = getSignalPlanarSensor(xtrack,Qin,D,S,N0,isNtype,isOxygenated,Vbias_2,ifluence,T, t,dt, N);
%             CCE_1(iter)=abs(min(Q_tot_1)/Qin)*100;
%             CCE_2(iter)=abs(min(Q_tot_2)/Qin)*100;
%             iter = iter + 1;
%             ifluence
%         end
%         
%         plot(fluence, CCE_1, 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '-');
%         hold on;
%         plot(fluence, CCE_2, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '-');
%         hold off;
%         title('Charge collection efficiency', 'FontWeight','bold','FontSize', 10);
%         xlabel('fluence [10^{12} N_{eq}/cm^2]', 'FontWeight','bold');
%         ylabel('CCE [%]', 'FontWeight','bold');
%         legend('V_{bias} = 1.1*V_{dep}', 'V_{bias} = V_{dep}+20V','Location', 'northeast');
%         grid on;
%         set(gca, 'GridLineStyle', '-');
%         ylim([0 100]);
%         hold off;
%         export_fig('CCE_irradiation.pdf');

%   plot signal as a function of time
%     [Q_tot Q_ind_e_vec Q_ind_h_vec y_e_vec y_h_vec v_e_vec v_h_vec E_q_y_e_vec E_q_y_h_vec E_w_y_e_vec E_w_y_h_vec Phi_w_y_e_vec Phi_w_y_h_vec] = getSignalPlanarSensor(xtrack,Qin,D,S,N0,isNtype,isOxygenated,Vbias_1,fluence,T, t,dt, N);
    
    [Q_tot Q_ind_e_vec Q_ind_h_vec y_e_vec y_h_vec v_e_vec v_h_vec E_q_y_e_vec E_q_y_h_vec E_w_y_e_vec E_w_y_h_vec Phi_w_y_e_vec Phi_w_y_h_vec] = getSignalPlanarSensor(xtrack,Qin,D,S,N0,isNtype,isOxygenated,Vbias_1,fluence,T, t,dt, N);

% plot induced charge as a function of time for all e-h pairs
    
%         plot(t, sum(Q_ind_e_vec,1)./1e5+sum(Q_ind_h_vec,1)./1e5, 'COLOR', 'black','LineWidth', 2, 'LineStyle', '-');
%         hold on;
%         plot(t, sum(Q_ind_e_vec,1)./1e5, 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '-');
%     	plot(t, sum(Q_ind_h_vec,1)./1e5, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '-');
%         set(gcf, 'Color', [1 1 1]);
%         title_str = sprintf('Induced charge,particle track,n-oxyg. silicon sensor,300 um,N_{eff}=%1.0f*10^{12}/cm³,V_{bias}=1.5 V_{dep}=%1.0f V', Neff, Vbias_1);
%         title(title_str, 'FontWeight','bold','FontSize', 10);
%         xlabel('time [ns]', 'FontWeight','bold');
%         ylabel('total induced charge [10 ke]', 'FontWeight','bold');
%         legend('e+h: Q_{ind}', 'e: Q_{ind}', 'h: Q_{ind}','Location', 'southwest');
%         grid on;
%         set(gca, 'GridLineStyle', '-');
%         ylim([-1 0]);
%         hold off;
%         export_fig('signal_track_N50.pdf');
    % plot the silicon detector with animated charge
%         deltaT = floor(size(t,2)/100);  %time step width for movie
%         frameIndex = 1;
%         iterations = 1:deltaT:floor(size(t,2));
%         plotPlanarDetector(D, S, 100, 10, y_e_vec(:,1), y_h_vec(:,1), deltaT./1000);
%         f = getframe(gcf);
%         [im,map] = rgb2ind(f.cdata,256,'nodither');
%         im(1,1,1,frameIndex) = 0;
%         for (j=iterations)
%             j
%             plotPlanarDetector(D, S, 100, 10, y_e_vec(:,j), y_h_vec(:,j), j./1000);
%             f = getframe(gcf);
%             im(:,:,1,frameIndex) = rgb2ind(f.cdata,map,'nodither');
%             M(frameIndex) = getframe;
%             frameIndex=frameIndex+1;
%         end
%         imwrite(im, map,'planar_3_pos.gif','DelayTime',0.005,'LoopCount',inf); %g443800´
    
%     % plot values as a function of time for one e-h pair
        plot(t, Q_ind_e_vec(1,:)./1e5+Q_ind_h_vec(1,:)./1e5, 'COLOR', 'black','LineWidth', 2, 'LineStyle', '-');
        hold on;
        plot(t, Q_ind_e_vec(2,:)./1e5+Q_ind_h_vec(2,:)./1e5, 'COLOR', 'black','LineWidth', 2, 'LineStyle', '--');
        plot(t, Q_ind_e_vec(3,:)./1e5+Q_ind_h_vec(3,:)./1e5, 'COLOR', 'black','LineWidth', 2, 'LineStyle', ':');
%         plot(t, y_e_vec(1,:)./100, 'COLOR', 'black','LineWidth', 1.2, 'LineStyle', '-');
%         plot(t, y_h_vec(1,:)./100, 'COLOR', 'black','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, v_e_vec(1,:)./100, 'COLOR', 'red','LineWidth', 1.2);
%         plot(t, v_h_vec(1,:)./100, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '--');
        plot(t, Q_ind_e_vec(2,:)./1e5, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '-');
        plot(t, Q_ind_h_vec(2,:)./1e5, 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '-');
%         plot(t, Q_ind_e_vec(2,:)./1e5, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, Q_ind_h_vec(2,:)./1e5, 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, E_q_y_e_vec(1,:), 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '-');
%         plot(t, E_q_y_h_vec(1,:), 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, E_w_y_e_vec(2,:).*100, 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '-');
%         plot(t, E_w_y_h_vec(2,:).*100, 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, Phi_w_y_e_vec(1,:), 'COLOR', 'green','LineWidth', 1.2, 'LineStyle', '-');
%         plot(t, Phi_w_y_h_vec(1,:), 'COLOR', 'green','LineWidth', 1.2, 'LineStyle', '--');

        ylim([-1 0]);
        xlim([0 5]);
        set(gcf, 'Color', [1 1 1]);
        title('Q_{meas} for Q_{ind} at different positions, 300 um, 1e12/cm^3, V_{bias}=1.5 V_{dep}', 'FontWeight','bold','FontSize', 10);
        xlabel('time [ns]', 'FontWeight','bold');
        ylabel('Q_{meas} [|Q_{in}|]', 'FontWeight','bold');
%         ylabel('differend quantities, see legend', 'FontWeight','bold');
        legend('1: Q_{in} at p electrode','2: Q_{in} middle of detector','3: Q_{in} at n electrode','electron signal for 2','hole signal for 2','Location', 'northeast');
%         legend('Q_{ind} [10 ke]','e: position [100 um]','h: position [100 um]','e: velocity [100 um/ns]','h: velocity [100 um/ns]', 'e: Q_{ind} [10 ke]','h: Q_{ind} [10 ke]', 'e: e-field [V/um]','h: e-field [V/um]' , 'e: weighting-pot.','h: weighting-pot.','Location', 'northeast');
        grid on;
        set(gca, 'GridLineStyle', '-');
        hold off;
        export_fig('signal_diff_pos.pdf');



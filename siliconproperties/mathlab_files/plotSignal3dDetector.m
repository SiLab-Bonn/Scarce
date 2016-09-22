%constants
    clear;
    T = 248;    %temperature [K]
    fluence = 0e4;    %radiation dose [10^12 Neq/cm^3]
    N0 = 1.00; %doping concentration [10^12 /cm^3]
    isNtype = 1; %type of the sensor
    isOxygenated = 1; %oxygenated sensor
    D = 200; % detector width [um]
    S = 70; %electrode distance [um]
    R = 5; %column radius [um]
    Qin = 1e5;  %deposited charge in electron hole pairs
%     xtrack = [0 0 0];
%     ytrack = [D/2 D/3 D/2];
    xtrack = linspace(-S/2+R,S/2-R,10);  %offset from the middle of the sensor pixel (x-direction) [um]
    xtrack = repmat(xtrack,[size(xtrack,2) 1]);
    xtrack = xtrack(:)';
    ytrack = linspace(D/2-20,D/2+20,floor(sqrt(size(xtrack,2))));  %offset from the bottom of the sensor pixel (y-direction) [um]
    ytrack = repmat(ytrack,1,size(ytrack,2));
    Neff = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated);
    Vbias_1 = -1.1*getDepletionVoltage(Neff, D)-20; %bias voltage [V];
    dt = 1e-4;   %time step of simulation [ns]
    t = 0:dt:1.0; %time points of simulation [ns]
%     t = 0;
    Resolution = 200;
        
%   plot signal as a function of time
    [Q_tot Q_ind_e_vec Q_ind_h_vec x_e_vec x_h_vec v_x_e_vec v_x_h_vec y_e_vec y_h_vec v_y_e_vec v_y_h_vec E_q_x_e_vec E_q_x_h_vec E_w_x_e_vec E_w_x_h_vec E_q_y_e_vec E_q_y_h_vec E_w_y_e_vec E_w_y_h_vec Phi_w_e_vec Phi_w_h_vec] = getSignal3dSensor(xtrack,ytrack,Qin,D,S,R,N0,isNtype,isOxygenated,Vbias_1,fluence,T, t,dt);
    
%     v_x_e_vec 
%     v_y_e_vec 
%         plot3dDetector(D, S, R, D/2, Resolution, 21, x_e_vec(:,1), y_e_vec(:,1), x_h_vec(:,1), y_h_vec(:,1), 0);

%     % plot the silicon detector with animated charge
        deltaT = floor(size(t,2)/100);  %time step width for movie
        frameIndex = 1;
        iterations = 1:deltaT:floor(size(t,2))

        plot3dDetector(D, S, R, D/2, Resolution, 21, x_e_vec(:,1), y_e_vec(:,1), x_h_vec(:,1), y_h_vec(:,1), deltaT.*dt);
        f = getframe(gcf);
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,frameIndex) = 0;
        for (j=iterations)
            j
            plot3dDetector(D, S, R, D/2, Resolution, 21, x_e_vec(:,j), y_e_vec(:,j), x_h_vec(:,j), y_h_vec(:,j), j.*dt);
            f = getframe(gcf);
            im(:,:,1,frameIndex) = rgb2ind(f.cdata,map,'nodither');
            M(frameIndex) = getframe;
            frameIndex=frameIndex+1;
        end
        imwrite(im, map,'3d_animation.gif','DelayTime',0.005,'LoopCount',inf); %g443800´
%     
%     % plot values as a function of time for one track
%         plot(t, Q_tot./1e5, 'COLOR', 'blue','LineWidth', 2, 'LineStyle', '-');   
% %         plot(t, Q_ind_e_vec(1,:)./1e5+Q_ind_h_vec(1,:)./1e5, 'COLOR', 'blue','LineWidth', 2, 'LineStyle', '-');   
%         hold on;
%         plot(t, x_e_vec(1,:)./100, 'COLOR', 'black','LineWidth', 1.2, 'LineStyle', '-');
%         plot(t, y_e_vec(1,:)./100, 'COLOR', 'black','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, x_h_vec(1,:)./100, 'COLOR', 'black','LineWidth', 1.2, 'LineStyle', '-.');
%         plot(t, y_h_vec(1,:)./100, 'COLOR', 'black','LineWidth', 1.2, 'LineStyle', ':');
%         plot(t, v_x_e_vec(1,:)./100, 'COLOR', 'red','LineWidth', 1.2);
%         plot(t, v_y_e_vec(1,:)./100, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, v_x_h_vec(1,:)./100, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '-.');
%         plot(t, v_y_h_vec(1,:)./100, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', ':');
% %         plot(t, sqrt((v_x_e_vec(1,:)).^2+(v_y_e_vec(1,:).^2))./100, 'COLOR', 'red','LineWidth', 1.2, 'LineStyle', '-.');
%         plot(t, E_w_x_e_vec(1,:).*100, 'COLOR', 'green','LineWidth', 1.2, 'LineStyle', '-');
%         plot(t, E_w_y_e_vec(1,:).*100, 'COLOR', 'green','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, E_w_x_h_vec(1,:).*100, 'COLOR', 'green','LineWidth', 1.2, 'LineStyle', '-.');
%         plot(t, E_w_y_h_vec(1,:).*100, 'COLOR', 'green','LineWidth', 1.2, 'LineStyle', '--');
%         plot(t, Phi_w_e_vec(1,:), 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '-.');
%         plot(t, Phi_w_h_vec(1,:), 'COLOR', 'blue','LineWidth', 1.2, 'LineStyle', '-.');
%         ylim([-1.5 2.5]);
%         xlim([0 3]);
%         set(gcf, 'Color', [1 1 1]);
%         title_str = sprintf('3d-silicon sensor, %1.0f um electrode distance, %1.0f um radius, charge at (x,y) = (%1.0f,%1.0f) um', S, R, xtrack(1), ytrack(1));
%         title(title_str, 'FontWeight','bold','FontSize', 10);
%         xlabel('time [ns]', 'FontWeight','bold');
%         ylabel('differend quantities, see legend', 'FontWeight','bold');
%         legend('Q_{ind} [10 ke]', 'x_{e}[100 um]','y_{e}[100 um]','x_{h}[100 um]','y_{h}[100 um]','vx_{e}[100 um/ns]','vy_{e}[100 um/ns]','vx_{h}[100 um/ns]','vy_{h}[100 um/ns]', 'e: w-field x [V/um]', 'e: w-field y [V/um]','h: w-field x [V/um]','h: w-field y [V/um]','e: w-pot [a.u.]','h: w-field [a.u.]','Location', 'northeast');
%         grid on;
%         set(gca, 'GridLineStyle', '-');
%         hold off;
% %         export_fig('signal_3d.pdf');

    % plot values as a function of time for even distributed tracks
        plot(t, mean(Q_ind_e_vec(:,:)+Q_ind_h_vec(:,:),1)./Qin, 'COLOR', 'red','LineWidth', 2, 'LineStyle', '-');
        hold on;
        plot(t, Q_ind_e_vec(1,:)./Qin+Q_ind_h_vec(1,:)./Qin, 'COLOR', 'red','LineWidth', 2, 'LineStyle', '--');
        plot(t, Q_ind_e_vec(2,:)./Qin+Q_ind_h_vec(2,:)./Qin, 'COLOR', 'red','LineWidth', 2, 'LineStyle', '-.');
        hold off;
        ylim([-1.0 0]);
        xlim([0 2]);
        set(gcf, 'Color', [1 1 1]);
        title_str = sprintf('3d-silicon sensor, %1.0f um electrode distance, %1.0f um radius, charges at (x,y) = (-%1.0f:%1.0f,-%1.0f:%1.0f) um', S, R, xtrack(1), ytrack(1));
        title(title_str, 'FontWeight','bold','FontSize', 10);
        xlabel('time [ns]', 'FontWeight','bold');
        ylabel('differend quantities, see legend', 'FontWeight','bold');
        legend('Q_{ind} total [10 ke]','Location', 'northeast');
        grid on;
        set(gca, 'GridLineStyle', '-');


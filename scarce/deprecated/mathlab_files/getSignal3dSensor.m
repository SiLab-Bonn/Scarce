% function to return the signal of a 3d silicon sensor as a function of
% the time. The parameters are:
%   xtrack: offset from the electrode mean [um]
%   Qin: deposited charge [e-h pairs]
%   D: sensor thickness [um]
%   S: electrode distance [um]
%   R: electrode radius [um]
%   N0: doping concentration without fluence [10^12 Neq/cm^3]
%   isNtype: = 1 for ntype sensor
%   isOxygenated: = 1 for oxygenated sensor
%   Vbias: bias of the sensor [V]
%   fluence: radiation dose [10^12 Neq/cm^3]
%   T: temerature [K]
%   t: time points of simulation [ns]
%   dt: time step of simulation [ns] CHOOSE SMALL ENOUGH!
%   N: number of quasi particles along the track

function [Q_tot Q_ind_e_vec Q_ind_h_vec x_e_vec x_h_vec v_x_e_vec v_x_h_vec y_e_vec y_h_vec v_y_e_vec v_y_h_vec E_q_x_e_vec E_q_x_h_vec E_w_x_e_vec E_w_x_h_vec E_q_y_e_vec E_q_y_h_vec E_w_y_e_vec E_w_y_h_vec Phi_w_e_vec Phi_w_h_vec] = getSignal3dSensor(xtrack,ytrack,Qin,D,S,R,N0,isNtype,isOxygenated,Vbias,fluence,T, t,dt)
    % constants
    e_h = Qin;  %deposited charge in electron hole pairs
    Neff = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated);
    tr_e = getTrappingTime(fluence,1); %trapping time of electrons [ns]
    tr_h = getTrappingTime(fluence,0); %trapping time of holes [ns]
    
    E = 1e6; % maximum electric field after strong irradiation [V/cm]
    v_max = getMobility(E, T, 1).*E; % maximum saturation velocity [cm/s]
 
    x_e = xtrack;   %positions of the electrons in x
    y_e = ytrack;   %positions of the electrons in y
    
    x_h = x_e; %positions of the electron holes in x    
    y_h = y_e; %positions of the electron holes in y    
    
    Q_ind_e = linspace(0,0,size(y_e,2));%start induced charge from electrons
    Q_ind_h = linspace(0,0,size(y_h,2));%start induced charge from holes
    index = 1;  %counting index to create plot at each time point index

    Q_ind_e_vec = zeros(size(y_e, 2), size(t, 2));
    Q_ind_h_vec = zeros(size(y_h, 2), size(t, 2));
    y_e_vec = zeros(size(y_e, 2), size(t, 2));
    y_h_vec = zeros(size(y_h, 2), size(t, 2));
    v_y_e_vec = zeros(size(y_e, 2), size(t, 2));
    v_y_h_vec = zeros(size(y_h, 2), size(t, 2));
    x_e_vec = zeros(size(y_e, 2), size(t, 2));
    x_h_vec = zeros(size(y_h, 2), size(t, 2));
    v_x_e_vec = zeros(size(y_e, 2), size(t, 2));
    v_x_h_vec = zeros(size(y_h, 2), size(t, 2));
    E_q_y_e_vec = zeros(size(y_e, 2), size(t, 2));
    E_q_y_h_vec = zeros(size(y_h, 2), size(t, 2));
    E_w_y_e_vec = zeros(size(y_e, 2), size(t, 2));
    E_w_y_h_vec = zeros(size(y_h, 2), size(t, 2));
    E_q_x_e_vec = zeros(size(x_e, 2), size(t, 2));
    E_q_x_h_vec = zeros(size(x_h, 2), size(t, 2));
    E_w_x_e_vec = zeros(size(x_e, 2), size(t, 2));
    E_w_x_h_vec = zeros(size(x_h, 2), size(t, 2));
    Phi_w_e_vec = zeros(size(y_e, 2), size(t, 2));
    Phi_w_h_vec = zeros(size(y_h, 2), size(t, 2));


    for(i = t)  %time loop
        % electric, weighting field a charge carrier position
        % electrons 
            [E_q_x_e E_q_y_e] = getWeightingField3d(x_e, y_e-D/2, S, R);   % electric field [V/um]
            Phi_w_e = getWeightingPotential3d(x_e, y_e-D/2, S, R);
            [E_w_x_e E_w_y_e] = getWeightingField3d(x_e, y_e-D/2, S, R);
        % holes
            [E_q_x_h E_q_y_h] = getWeightingField3d(x_h, y_h-D/2, S, R);   % electric field [V/um]
            Phi_w_h = getWeightingPotential3d(x_h, y_h-D/2, S, R);
            [E_w_x_h E_w_y_h] = getWeightingField3d(x_h, y_h-D/2, S, R);

        % movement
            % electrons
%             v_e_x = -E_q_x_e.*getMobility(E_q_x_e.*1e5, T, 1); % velocity in x [um/ns]
%             v_e_y = -E_q_y_e.*getMobility(E_q_y_e.*1e5, T, 1); % velocity in y [um/ns]
            % simple assumption that the velocity is constant
            E_q_norm = sqrt(E_q_x_e.^2+E_q_y_e.^2);
            v_e_x = -E_q_x_e./E_q_norm.*v_max.*1e-5; %um/ns
            v_e_y = -E_q_y_e./E_q_norm.*v_max.*1e-5; %um/ns

            dx_e = v_e_x .* dt;
            dy_e = v_e_y .* dt;
            x_e = x_e + dx_e;
            y_e = y_e + dy_e;
            
            % set limits
            x_e(isnan(y_e))=NaN;
            y_e(isnan(x_e))=NaN;
            E_w_x_e(isnan(x_e)) = NaN;
            E_w_y_e(isnan(x_e)) = NaN;
            E_q_x_e(isnan(x_e)) = NaN;
            E_q_y_e(isnan(x_e)) = NaN;
            v_e_x(isnan(x_e)) = NaN;
            v_e_y(isnan(x_e)) = NaN; 
            
            % holes
%             v_h_x = -E_q_x_h.*getMobility(E_q_x_h.*1e5, T, 0); % velocity in x [um/ns]
%             v_h_y = -E_q_y_h.*getMobility(E_q_y_h.*1e5, T, 0); % velocity in y [um/ns]
            E_q_norm = sqrt(E_q_x_h.^2+E_q_y_h.^2);
            v_h_x = E_q_x_h./E_q_norm.*v_max.*1e-5; %um/ns
            v_h_y = E_q_y_h./E_q_norm.*v_max.*1e-5; %um/ns
            
            dx_h = v_h_x .* dt;
            dy_h = v_h_y .* dt;
            x_h = x_h + dx_h;
            y_h = y_h + dy_h;
            
            % set limits
            x_h(isnan(y_h))=NaN;
            y_h(isnan(x_h))=NaN;
            E_w_x_h(isnan(x_h)) = NaN;
            E_w_y_h(isnan(x_h)) = NaN;
            E_q_x_h(isnan(x_h)) = NaN;
            E_q_y_h(isnan(x_h)) = NaN;
            v_h_x(isnan(x_h)) = NaN;
            v_h_y(isnan(x_h)) = NaN; 

        % induced charge calculation
            % electrons
            dQ_e = linspace(0,0,size(y_e, 2));
            dQ_e(~isnan(x_e)) = e_h.*exp(-i./tr_e).*(E_w_x_e(~isnan(x_e)).*dx_e(~isnan(x_e))+E_w_y_e(~isnan(x_e)).*dy_e(~isnan(x_e)));

            % holes
            dQ_h = linspace(0,0,size(y_h, 2));
            dQ_h(~isnan(x_h)) = -e_h.*exp(-i./tr_h).*(E_w_x_h(~isnan(x_h)).*dx_h(~isnan(x_h))+E_w_y_h(~isnan(x_h)).*dy_h(~isnan(x_h)));

            Q_ind_e = Q_ind_e + dQ_e;
            Q_ind_h = Q_ind_h + dQ_h;

        % data for plotting
        Q_ind_e_vec(:,index) = Q_ind_e;
        Q_ind_h_vec(:,index) = Q_ind_h;
        x_e_vec(:,index) = x_e;
        v_x_e_vec(:,index) = v_e_x;
        x_h_vec(:,index) = x_h;
        v_x_h_vec(:,index) = v_h_x;
        y_e_vec(:,index) = y_e;
        v_y_e_vec(:,index) = v_e_y;
        y_h_vec(:,index) = y_h;
        v_y_h_vec(:,index) = v_h_y;
        E_q_x_e_vec(:,index) = E_q_x_e;
        E_q_x_h_vec(:,index) = E_q_x_h;
        E_w_x_e_vec(:,index) = E_w_x_e;
        E_w_x_h_vec(:,index) = E_w_x_h;       
        E_q_y_e_vec(:,index) = E_q_y_e;
        E_q_y_h_vec(:,index) = E_q_y_h;
        E_w_y_e_vec(:,index) = E_w_y_e;
        E_w_y_h_vec(:,index) = E_w_y_h;
        Phi_w_e_vec(:,index) = Phi_w_e;
        Phi_w_h_vec(:,index) = Phi_w_h;

        index = index + 1;  %increase plot data index
    end;

    Q_tot = sum(Q_ind_e_vec,1)+sum(Q_ind_h_vec,1);   
end %getSignal3dSensor




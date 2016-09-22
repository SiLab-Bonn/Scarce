% function to return the signal of a planar silicon sensor as a function of
% the time. CONSTANT VELOCITY ASSUMED. The parameters are:
%   xtrack: offset from the electrode mean [um]
%   Qin: deposited charge [e-h pairs]
%   D: sensor thickness [um]
%   S: electrode width [um]
%   N0: doping concentration without fluence [10^12 Neq/cm^3]
%   isNtype: = 1 for ntype sensor
%   isOxygenated: = 1 for oxygenated sensor
%   Vbias: bias of the sensor [V]
%   fluence: radiation dose [10^12 Neq/cm^3]
%   T: temerature [K]
%   t: time points of simulation [ns]
%   dt: time step of simulation [ns] CHOOSE SMALL ENOUGH!
%   N: number of quasi particles along the track

function [Q_tot Q_ind_e_vec Q_ind_h_vec y_e_vec y_h_vec v_e_vec v_h_vec E_q_y_e_vec E_q_y_h_vec E_w_y_e_vec E_w_y_h_vec Phi_w_y_e_vec Phi_w_y_h_vec] = getSignalPlanarSensorSimple(xtrack,Qin,D,S,N0,isNtype,isOxygenated,Vbias,fluence,T, t,dt, N)
    % constants
    Q_e = 1.602e-19;   %electron charge [C]
    e_h = Qin;  %deposited charge in electron hole pairs
    x_0 = xtrack;  %offset from the middle of the sensor pixel (x-direction) [um]
    Neff = getEffAcceptorConentration(fluence, N0, isNtype, isOxygenated);
    tr_e = getTrappingTime(fluence,1); %trapping time of electrons [ns]
    tr_h = getTrappingTime(fluence,0); %trapping time of holes [ns]
    
    E = 1e7; % maximum electric field after strong irradiation [V/cm]
    v_max = getMobility(E, T, 1).*E; % maximum saturation velocity [cm/s]

    e_h = e_h/N;

    %start positions of the electrons
    if(N==1)
        y_e = D/2;
    else
        y_e = linspace(0,D,N);
    end
    y_h = y_e;  %start positions of the electron holes
    x_e = linspace(x_0,x_0,size(y_e,2));%positions of the electrons in x
    x_h = linspace(x_0,x_0,size(y_e,2));%positions of the electron holes in x
    Q_ind_e = linspace(0,0,size(y_e,2));%start induced charge from electrons
    Q_ind_h = linspace(0,0,size(y_h,2));%start induced charge from holes
    index = 1;  %counting index to create plot at each time point index

    Q_ind_e_vec = zeros(size(y_e, 2), size(t, 2));
    Q_ind_h_vec = zeros(size(y_h, 2), size(t, 2));
    y_e_vec = zeros(size(y_e, 2), size(t, 2));
    y_h_vec = zeros(size(y_h, 2), size(t, 2));
    v_e_vec = zeros(size(y_e, 2), size(t, 2));
    v_h_vec = zeros(size(y_h, 2), size(t, 2));
    E_q_y_e_vec = zeros(size(y_e, 2), size(t, 2));
    E_q_y_h_vec = zeros(size(y_h, 2), size(t, 2));
    E_w_y_e_vec = zeros(size(y_e, 2), size(t, 2));
    E_w_y_h_vec = zeros(size(y_h, 2), size(t, 2));
    Phi_w_y_e_vec = zeros(size(y_e, 2), size(t, 2));
    Phi_w_y_h_vec = zeros(size(y_h, 2), size(t, 2));

    for(i = t)  %time loop
        % electric, weighting field a charge carrier position
        % electrons 
            E_q_y_e = getElectricFieldPlanar(y_e, Vbias, Neff, D);    % electric field [V/um]
            Phi_w_e = getWeightingPotentialPlanar(x_e, y_e, D, S);
            [E_w_x_e E_w_y_e] = getWeightingFieldPlanar(x_e, y_e, D, S);
        % holes
            E_q_y_h = getElectricFieldPlanar(y_h, Vbias, Neff, D);    % electric field [V/um]
            Phi_w_h = getWeightingPotentialPlanar(x_h, y_h, D, S);
            [E_w_x_h E_w_y_h] = getWeightingFieldPlanar(x_h, y_h, D, S);

        % movement
            % electrons
%             v_e = -E_q_y_e.*getMobility(E_q_y_e.*1e5, T, 1); % velocity [um/ns]
            v_e = linspace(1,1,size(E_q_y_e,2)).*v_max.*1e-5;% velocity [um/ns]
            dy_e = v_e .* dt;
            y_e = y_e + dy_e;
            % set limits
            E_w_y_e(y_e>D) = 0;
            E_q_y_e(y_e>D) = 0;
            v_e(y_e>=D) = 0;
            y_e(y_e>D) = D;      

            % holes
%             v_h = E_q_y_h.*getMobility(E_q_y_h.*1e5, T, 0); % velocity [um/ns]
            v_h = -linspace(1,1,size(E_q_y_h,2)).*v_max.*1e-5;
            dy_h = v_h .* dt;
            y_h = y_h + dy_h;
            % set limits
            E_w_y_h(y_h<0) = 0;
            E_q_y_h(y_h<0) = 0;
            v_h(y_h<=0) = 0;
            y_h(y_h<0) = 0;

        % induced charge calculation
            % electrons
            dQ_e = linspace(0,0,size(y_e, 2));
            dQ_e(y_e<D) = e_h.*exp(-i./tr_e).*E_w_y_e(y_e<D).*dy_e(y_e<D);

            % holes
            dQ_h = linspace(0,0,size(y_h, 2));
            dQ_h(y_h>0) = -e_h.*exp(-i./tr_h).*E_w_y_h(y_h>0).*dy_h(y_h>0);

            Q_ind_e = Q_ind_e + dQ_e;
            Q_ind_h = Q_ind_h + dQ_h;

        % data for plotting
        Q_ind_e_vec(:,index) = Q_ind_e;
        Q_ind_h_vec(:,index) = Q_ind_h;
        y_e_vec(:,index) = y_e;
        v_e_vec(:,index) = v_e;
        y_h_vec(:,index) = y_h;
        v_h_vec(:,index) = v_h;
        E_q_y_e_vec(:,index) = E_q_y_e;
        E_q_y_h_vec(:,index) = E_q_y_h;
        E_w_y_e_vec(:,index) = E_w_y_e;
        E_w_y_h_vec(:,index) = E_w_y_h;
        Phi_w_y_e_vec(:,index) = Phi_w_e;
        Phi_w_y_h_vec(:,index) = Phi_w_h;

        index = index + 1;  %increase plot data index
    end;

    Q_tot = sum(Q_ind_e_vec,1)+sum(Q_ind_h_vec,1);   
end %getSignalPlanarSensor




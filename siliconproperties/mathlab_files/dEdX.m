clear;
%constants
    z = 1;                      %charge of incident particle [e]
    mp = 938.27231;             %mass of a proton [MeV/c^2]
    re = 2.817940325*10^-13;	%electron radius [cm]
    Z = 14;                     %atomic number [#nucleons]
    A = 28.0855;                %atomic mass of silicon absorber [amu = g/mol]
    Na = 6.0221415*10^23;       %avogrado number
    me = 0.510998918;           %mass of an electron [MeV/c^2]
    I = 173*10^(-6);            %mean excitation energy [MeV]
    K = 4*pi*Na*re^2*me;        %constant(?) [MeV]
    Alpha = 1/137.03599911;     %fine structure constant
    
    %density correction parameters for silicon, Atomic data and nuclear
    %tables 30, 261-271 (1984)
    Omega_p = 31.06;            %plasma energy in silicon [eV]
    C_d = -4.4351;              %density correction parameter
    A_d = 0.14921;              %density correction parameter
    M_d = 3.2546;               %density correction parameter
    X_0_d = 0.2014;             %density correction parameter
    X_1_d = 2.8715;             %density correction parameter
    
    %restricted Eloss, thin absorbers correction parameters
    T_cut = 19e-3;              %[MeV] for 200 um silicon from Nucl. Phys. B288 (1987) 681-716
    
    %MP Eloss, Landau-Vavilov Bichsel
    detector_density = 2.3290;  %density of detector material [g cm^-3]
    det_width = 230;            %detector width [um]
    mat_length = det_width*1e-4*2.3290; %material length [g/cm^2]
    j = 0.198;                    %parameter [no unit]
    
    %MP electron hole pairs measured
    electron_hole_pair = 3.61;  %energy neede to create electron hole pair [eV]
    
%input
    betagamma = logspace(-1,3.3,1000);    %input
    %some derivations of the input to simplify the formulars
    betagamma2 = betagamma.*betagamma;  %square betagamma
    gamma = (betagamma2+1).^0.5;        %calculate gamma from betagamma squared
    beta = (1-1./(1+betagamma2)).^0.5;  %calculate beta from betagamma squared           
    beta2 = beta.*beta;                 %square beta
    %beta = (1-1./(gamma.*gamma)).^0.5; %cross check

%calculate dEdx
    %maximum energy transfer possible in single collision
    Tmax = (2.0*me.*betagamma2)./(1.0 + 2.0.*gamma.*(me./mp) + (me./mp).^2);
    
    %density effect correction(PDG 30. Passage of particles through
    %matter)
    X = log10(betagamma);   %derivation of the input to simplify the formular
    %calculate delta with three cases
    delta=zeros(1, size(betagamma,2));                  %else case: zero for X < X_0
    delta(X>=X_1_d) = 2.*log(10).*X(X>=X_1_d) + C_d;    %case X >= X1: 2*ln 10 * X - C
    delta(X>=X_0_d & X<X_1_d) = 2.*log(10).*X(X>=X_0_d & X<X_1_d) + C_d + A_d.*(X_1_d - X(X>=X_0_d & X<X_1_d)).^M_d;                %case X>=X_0_d & X<X_1_d: 2*log(10)*X + C_d + A_d*(X_1_d - X)^M_d
    
    %MP Eloss correction
    kappa = K./2.*Z./A.*mat_length./beta2;
    
    %semiempiric Bethe Bloch formular
    dedx = (-K*z.^2*Z)./(A.*beta2).*(0.5.*log((2.*me.*betagamma2.*Tmax)./I.^2)-beta2);                                              %dEdx
    dedx_density = (-K*z.^2*Z)./(A.*beta2).*(0.5.*log((2.*me.*betagamma2.*Tmax)./I.^2)-beta2-delta./2);                             %dEdx with density effect (with silicon parameters)
    dedx_density_rest = (-K*z.^2*Z)./(A.*beta2).*(0.5.*log((2.*me.*betagamma2.*T_cut)./I.^2)-beta2./2.*(1+T_cut./Tmax)-delta./2);   %dEdx with density effect and restriction on measured dEdx in thin absorbers
    dedx_density_MPV = kappa.*(log(2.*me.*betagamma2./I)+log(kappa./I)+j-beta2-delta)./(mat_length);

%MP charge deposit
    dedx_Q = dedx_density_MPV./electron_hole_pair.*1e3.*mat_length; %factor 1e3 for scaling to ke

%plot dEdX
    h = semilogx(betagamma,-dedx, 'LineWidth', 2);
    hold on;
    semilogx(betagamma,-dedx_density, 'LineWidth', 2, 'color', 'red');
    semilogx(betagamma,-dedx_density_rest, 'LineWidth', 2, 'color', 'green');
    semilogx(betagamma,dedx_density_MPV, 'LineWidth', 2, 'color', 'black');
    %semilogx(betagamma,delta, 'LineWidth', 2, 'color', 'black');
    hold off;
    set(gca, 'YMinorTick','on');
    title('Average stopping power (energy loss per material distance) ', 'FontWeight','bold','FontSize', 10);
    xlabel('\beta\gamma = p/(Mc)', 'FontWeight','bold');
    ylabel('average stopping power [MeV g^{-1} cm^2]', 'FontWeight','bold');
    xlim([0.3 2000]);
    ylim([0 5]);
    set(gcf, 'Color', [1 1 1]);
    hleg1 = legend('Bethe Bloch','+ density correction','+ restricted Eloss','+ MP Eloss');
    grid on
%     xlabh = get(gca,'XLabel');
%     x2 = axes('Position',get(gca,'Position')-[0.0 0.05 0 0]);
%     set(x2,'ytick',[1 2 3 4 5], 'yticklabel',{})
    
    %axes('Position',[0 0 1 1]);
    %export_fig('dEdx.pdf');
    figure;
    semilogx(betagamma,dedx_Q, 'LineWidth', 2);
    min(dedx_Q)
    set(gca, 'YMinorTick','on');
    title_str = sprintf('Most propable charge deposited in %d um silicon layer', det_width);
    title(title_str, 'FontWeight','bold','FontSize', 10);
    xlabel('\beta\gamma = p/(Mc)', 'FontWeight','bold');
    ylabel('most propable charge deposited [ke]', 'FontWeight','bold');
    set(gca, 'YMinorTick','on');
    xlim([0.55 2000]);
    ylim([0 60]);
    grid on
    set(gcf, 'Color', [1 1 1]);
    export_fig('ChargeDeposited.pdf');
    min(dedx_Q)

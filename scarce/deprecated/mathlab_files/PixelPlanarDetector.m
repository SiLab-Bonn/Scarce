% Weighting field for pixel/strip planar detector calculations
    clear;

% constant
    DETWIDTH = 1;
    STRIPWIDTH = 0.6;   % in units of the DETWIDTH
    NPIXEL = 1; %numbers of pixelelectrodes to plot
    %DELTEPIXEL = [0 2 5 4];    %position of pixel in STRIPWIDTH to delete

% options
    RESOLUTION = 10;
    NFIELDPLOTARROWS = 20;
    SIZEPLOTARROWS = 0.5;
    XPLOTPART = 0.5;    %the part of the total grid size in x that is plotted

% calculation grid size
    GRIDSIZEX = DETWIDTH/RESOLUTION;
    GRIDSIZEY = DETWIDTH/RESOLUTION;
    x = -4*DETWIDTH:GRIDSIZEX:4*DETWIDTH;
    y = 0:GRIDSIZEY:DETWIDTH;

% variable transformation
    [xx,yy] = meshgrid(x,y);

% weighting potential and field of one pixel
    phi_w  = getWeightingPotentialPlanar(xx, yy, DETWIDTH, STRIPWIDTH);
    [Ex_w Ey_w] = gradient(-phi_w, GRIDSIZEX, GRIDSIZEY);
    
    %calculate total potential from all pixels (NPIXEL) by superposition
    phi_tot = zeros(size(y,2),size(x,2));
    Ex_tot = zeros(size(y,2),size(x,2));
    Ey_tot = zeros(size(y,2),size(x,2));
    Nshiftbins = floor(size(y,2)*STRIPWIDTH); %shift by the pixel electrode bins
    for i = 0:NPIXEL-1
        if mod(i,2) == 0    %add electrode on right side from existing
            phi_tot = phi_tot + circshift(phi_w, [0 (i-i/2)*Nshiftbins]);
            Ex_tot = Ex_tot + circshift(Ex_w, [0 (i-i/2)*Nshiftbins]);
            Ey_tot = Ey_tot + circshift(Ey_w, [0 (i-i/2)*Nshiftbins]);
        else    %add electrode on left side from existing
            phi_tot = phi_tot + circshift(phi_w, [0 -(i+1)/2*Nshiftbins]);
            Ex_tot = Ex_tot + circshift(Ex_w, [0 -(i+1)/2*Nshiftbins]);
            Ey_tot = Ey_tot + circshift(Ey_w, [0 -(i+1)/2*Nshiftbins]);
        end
    end
    Eabs_tot = sqrt(Ex_tot.^2 + Ey_tot.^2); % The length of each field vector

%plot results
    %figure; %create a new figure window
    %resize the results to the plot range
    phi_wquiver = phi_tot(:, floor((1-XPLOTPART)*size(phi_w,2)/2+1): floor((1+XPLOTPART)*size(phi_w,2)/2+1));
    Ex_wquiver = Ex_tot(:, floor((1-XPLOTPART)*size(Ex_tot,2)/2+1): floor((1+XPLOTPART)*size(Ex_tot,2)/2+1));
    Ey_wquiver = Ey_tot(:, floor((1-XPLOTPART)*size(Ey_tot,2)/2+1): floor((1+XPLOTPART)*size(Ey_tot,2)/2+1));
    Eabs_totquiver = Eabs_tot(:, floor((1-XPLOTPART)*size(Ey_tot,2)/2+1): floor((1+XPLOTPART)*size(Ey_tot,2)/2+1));
    yquiver = yy(:, floor((1-XPLOTPART)*size(yy,2)/2+1): floor((1+XPLOTPART)*size(yy,2)/2+1));
    xquiver = xx(:, floor((1-XPLOTPART)*size(xx,2)/2+1): floor((1+XPLOTPART)*size(xx,2)/2+1));
    phi_wquiver(phi_wquiver>1) = 1; %overwrite numerical instabilities
    %plot the potential
    imagesc([min(x)*XPLOTPART max(x)*XPLOTPART],[0 DETWIDTH], phi_wquiver);  %plot the weighting potential
    set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','normal','ydir','normal','dataAspectRatio',[1 1 1]);
    colorbar;
    hold on; % Add subsequent plots to the image
    xlabel('sensor width');
    ylabel('sensor heigth');
    title_str = sprintf('Weighting Potential and Weighting E-Field direction for a planar sensor with %d shorted pixel', NPIXEL);
    title(title_str, 'FontWeight','bold');
    if mod(i,2) == 0
        fill([-STRIPWIDTH*(NPIXEL)/2 -STRIPWIDTH*(NPIXEL)/2 STRIPWIDTH*(NPIXEL)/2 STRIPWIDTH*(NPIXEL)/2], [-0.005 0.005 0.005 -0.005], 'white'); %rectangle to show the pixel electrode
    else
        fill([-STRIPWIDTH*(NPIXEL+1)/2 -STRIPWIDTH*(NPIXEL+1)/2 STRIPWIDTH*(NPIXEL-1)/2 STRIPWIDTH*(NPIXEL-1)/2], [-0.005 0.005 0.005 -0.005], 'white'); %rectangle to show the pixel electrode
    end;
    %fill([-1 -1 1 1], [-0.05 0.05 0.05 -0.05], 'white'); %rectangle to show the pixel electrode
    fill([min(x) min(x) max(x) max(x)], [1 0.99 0.99 1], 'white'); %rectangle to show the frontside electrode
    %change the E field arrays to plot less arrows in the field plot
        arrowsarraydx = floor(size(xquiver,1)/NFIELDPLOTARROWS);
        arrowsarraydy = floor(size(yquiver,1)/NFIELDPLOTARROWS);
        Ex_wquiver = Ex_wquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        Ey_wquiver = Ey_wquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        xquiver = xquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        yquiver = yquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        Eabs_totquiver = Eabs_totquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        Ex_wquiver(Ex_wquiver>1) = NaN;
        Ex_wquiver(Ex_wquiver<-1) = NaN;
        Ey_wquiver(Ey_wquiver>1) = NaN;
        Ey_wquiver(Ey_wquiver<-1) = NaN;
    quiver(xquiver, yquiver, Ex_wquiver./Eabs_totquiver, Ey_wquiver./Eabs_totquiver,SIZEPLOTARROWS, 'Color','r');
    %quiver(xquiver, yquiver, Ex_wquiver, Ey_wquiver, 10);
    set(gcf, 'Position', [100 100 100+750 100+150], 'Color', [1 1 1]);
    hold off;           % Any subsequent plotting will overwrite the image!
    %export_fig '2dpixel5.pdf';
    %new windows for potential/E-Field at x=0
    figure;
    plot(y, phi_tot(:,floor(size(phi_tot,2)/2)+1), 'LineWidth',2, 'Color','blue');    %plot total potential at x = 0 
    hold on;
    %plot(-gradient(phi_w(:,floor(size(phi_w,2)/2)+1)));
    xlabel('y position [total sensor thickness]');
    ylabel('potential (blue), Electric field (red) [a.u.]');
    plot(y, Ey_tot(:,floor(size(Ey_tot,2)/2)+1), 'LineWidth',2, 'Color','red');      %plot Ey at x = 0
    plot(y, Ex_tot(:,floor(size(Ex_tot,2)/2)+1), 'LineWidth',2, 'Color','red');      %plot Ex at x = 0
    title_str = sprintf('Potential (blue) and E-Field in x,y (red) at x = 0, # pixel = %d', NPIXEL);
    title(title_str, 'FontWeight','bold');
    set(gcf, 'Color', [1 1 1]);
    hold off;
    %export_fig '1dpixel5.pdf';
    
   



function plotPlanarDetectorEfield(Vbias, D, S, R, N, y_e, y_h, time)
% The parameters are:
%   D: sensor thickness [um]
%   S: electrode width [um]
%   R: resolution (pixel per detector width) of the plot
%   N: = numbers of arrows fopr the E-Field plot
%   y_e, y_h: y position of the electron/holes
%   time: time of the simulation quoted in the title

Neff = 1; %1e12/cm³
% options
    
    SIZEPLOTARROWS = 0.7;

% calculation grid size
    GRIDSIZEX = D/R;
    GRIDSIZEY = D/R;
    x = -D:GRIDSIZEX:D;
    y = 0:GRIDSIZEY:D;

% variable transformation
    [xx,yy] = meshgrid(x,y);

% weighting potential and field of one pixel
    phi  = getWeightingPotentialPlanar(xx, yy, D*100, S);
    [Ex_w Ey_w] = getWeightingFieldPlanar(xx,yy, D,S);
    
    phi = getPotentialPlanar(yy, Vbias, Neff, D);
    Ey_w = -getElectricFieldPlanar(yy, Vbias, Neff, D);
    phi = phi-Vbias;
%     Ex_w = zeros(1,size(Ey_w,2));
    
%     size(Ex_w)
%     size(Ey_w)
    
%     phi_w = phi_w.*Vbias-Vbias;
%     phi_w
    E_abs_w = sqrt(Ex_w.^2 + Ey_w.^2); % The length of each field vector

%plot results
    %figure; %create a new figure window
    %resize the results to the plot range
    Ex_wquiver = Ex_w;
    Ey_wquiver = Ey_w;
    Eabs_totquiver = E_abs_w;
    yquiver = yy;
    xquiver = xx;
%     phi_w(phi_w>1) = 1; %overwrite numerical instabilities
    %plot the potential
    imagesc([min(x) max(x)],[0 D], phi);  %plot the weighting potential
    set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','normal','ydir','normal','dataAspectRatio',[1 1 1]);
    colorbar;
    colormap Gray;
    colormap(flipud(colormap));
    hold on; % Add subsequent plots to the image
    xlabel('sensor width');
    ylabel('sensor heigth'); 
    title_str = sprintf('Planar pixel sensor electric field, V_{bias} = -%1.0f, effective doping concentration = %1.0f*10^{12} cm^{-3}', Vbias, Neff);
    title(title_str, 'FontWeight','bold');
    
    %change the E field arrays to plot less arrows in the field plot
        arrowsarraydx = floor(size(xquiver,1)/N);
        arrowsarraydy = floor(size(yquiver,1)/N);
        Ex_wquiver = Ex_wquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        Ey_wquiver = Ey_wquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        xquiver = xquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        yquiver = yquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        Eabs_totquiver = Eabs_totquiver(1:arrowsarraydx:end, 1:arrowsarraydy:end);
        Ex_wquiver(Ex_wquiver>1) = NaN;
        Ex_wquiver(Ex_wquiver<-1) = NaN;
        Ey_wquiver(Ey_wquiver>1) = NaN;
        Ey_wquiver(Ey_wquiver<-1) = NaN;
    quiver(xquiver, yquiver, Ex_wquiver, Ey_wquiver,SIZEPLOTARROWS, 'Color','r');
    %quiver(xquiver, yquiver, Ex_wquiver, Ey_wquiver, 10);
    set(gcf, 'Position', [100 100 100+750 100+150], 'Color', [1 1 1]);
    %draw the charge carriers
        %holes
%         y_h(y_h<0) = 0;
%         plot(0,y_h,'r.','Color','red','MarkerSize',15);
%         %electrons
%         y_e(y_e>D) = D;
%         plot(0,y_e,'r.','Color','blue','MarkerSize',10);
    %draw electrodes
    fill([-S/2 -S/2 S/2 S/2], [1.05*D 0.98*D 0.98*D 1.05*D], 'red'); %rectangle to show the pixel electrode
    fill([min(x) min(x) max(x) max(x)], [-0.05*D 0.02*D 0.02*D -0.02*D], 'blue'); %rectangle to show the frontside electrode
    hold off;           % Any subsequent plotting will overwrite the image!
%      caxis([0 1]);
    export_fig 'planar_sensor_efield.pdf';
end
    
   



function plot3dDetectorEfield(Vbias, D, S, C, yOff, R, N, x_e, y_e, x_h, y_h, time)
% The parameters are:
%   D: sensor area in y plot [um]
%   S: electrode distance [um]
%   C: electrode radius [um]
%   yOff: electrode offset in y [um]
%   R: resolution (pixel per detector width) of the plot
%   N: = numbers of arrows for the E-Field plot
%   x_e, y_e, x_h, y_h: position of the electron/holes
%   time: time of the simulation quoted in the title
    
    SIZEPLOTARROWS = 0.7;

% calculation grid size
    GRIDSIZEX = S/R;
    GRIDSIZEY = D/R;
    x = -D:GRIDSIZEX:D;
    y = 0:GRIDSIZEY:D;
%     x
%     y
% variable transformation
    [xx,yy] = meshgrid(x,y);

% weighting potential and field of one pixel
    phi_w  = getWeightingPotential3d(xx, yy-yOff, S, C);
    [Ex_w Ey_w] = getWeightingField3d(xx,yy-yOff, S, C);
%     [Ex_w Ey_w] = getWeightingFieldPlanar(xx,yy-yOff,D,S);
    E_abs_w = sqrt(Ex_w.^2 + Ey_w.^2); % The length of each field vector

    phi_w = phi_w.*Vbias-Vbias;
    
%plot results
    %figure; %create a new figure window
    %resize the results to the plot range
    Ex_wquiver = Ex_w;
    Ey_wquiver = Ey_w;
    Eabs_totquiver = E_abs_w;
    yquiver = yy;
    xquiver = xx;
    phi_w(phi_w>Vbias) = Vbias; %overwrite numerical instabilities
    %plot the potential
    imagesc([min(x) max(x)],[0 D], phi_w);  %plot the weighting potential
    set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','normal','ydir','normal','dataAspectRatio',[1 1 1]);
    colorbar;
    colormap Gray;
    colormap(flipud(colormap));
    hold on; % Add subsequent plots to the image
    xlabel('sensor x');
    ylabel('sensor y'); 
    title_str = sprintf('3d pixel sensor with electric field, V_{bias} = -%1.0f V, no space charge', Vbias);
    title(title_str, 'FontWeight','bold');
    
    %draw cylinder
%     ang=0:0.01:2*pi; 
%     xp=C*cos(ang);
%     yp=C*sin(ang); 
%     plot(-S/2+xp,yOff+yp,'Color','red');
%     plot(S/2+xp,yOff+yp,'Color','blue');
    
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
%     %quiver(xquiver, yquiver, Ex_wquiver, Ey_wquiver, 10);
    set(gcf, 'Position', [100 100 100+750 100+150], 'Color', [1 1 1]);
    %draw the charge carriers
        %holes
%         plot(x_h,y_h,'r.','Color','red','MarkerSize',15);
        %electrons
%         plot(x_e,y_e,'r.','Color','blue','MarkerSize',10);
%     streamslice(xx, yy, Ex_w,Ey_w);
%     contourf(x,y,phi_w);
%      contour(x,y,phi_w, 'EdgeColor','none');
%     caxis([0 Vbias]);
    filledCircle([-S/2 yOff], C, 100, 'r');
    filledCircle([S/2 yOff], C, 100, 'b');
    daspect([1 1 1]);   %filled circle does that wrong
    hold off;           % Any subsequent plotting will overwrite the image!
    export_fig '3d_sensor_Efield.pdf';
end
    
   



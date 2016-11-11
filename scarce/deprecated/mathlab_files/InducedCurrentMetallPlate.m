function [Q_ind_tot] = InducedCurrentMetallPlate(X0)
    % Calculate/plot the induced charge into a grounded infinite metall plate
    % at x = 0 with z = 0 by using two mirror charges at x = X0 and -X0

    % constants
        X0 = X0; %+- distance of charges to the plate at t = t0 
        Q = -1;  %charge
        EPSILON_0 = 1;

    % options
        NFIELDPLOTARROWS = 50;
        SIZEPLOTARROWS = 0.5;

    % calculation grid size
        DIMENSIONSY = 10;
        DIMENSIONSZ = 10;
        GRIDSIZEY = 1000;
        GRIDSIZEZ = 1000;
        y = -DIMENSIONSY/2:DIMENSIONSY/GRIDSIZEY:DIMENSIONSY/2;
        z = -DIMENSIONSZ/2:DIMENSIONSZ/GRIDSIZEZ:DIMENSIONSZ/2;

    % variable transformation
        [yy,zz] = meshgrid(y,z);

    % potentials of the point charges
%         phi_1 = Q./(4*pi*EPSILON_0) .* 1./(sqrt(yy.^2+(xx+X0).^2));
%         phi_2 = -Q./(4*pi*EPSILON_0) .* 1./(sqrt(yy.^2+(xx-X0).^2));
%         phi_tot = phi_1 + phi_2;
%         %mesh(phi_tot);
%         [Ex_tot Ey_tot] = gradient(-phi_tot);
%         Eabs_tot = sqrt(Ex_tot.^2 + Ey_tot.^2); % The length of each field vector
        
    % charge surface density at metall plate (x=0)
        E_plate = Q.*X0./(2*pi*EPSILON_0*(X0.^2+yy.^2+zz.^2).^(3/2));
        
        Q_ind = EPSILON_0 .* E_plate;
        imagesc([-3 3], [-3 3], Q_ind, [-0.5 0]);  %plot the potential :, floor(GRIDSIZEX/2):GRIDSIZEX
        set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','normal','ydir','normal','dataAspectRatio',[1 1 1]);     
        colorbar;
        xlabel('y');
        ylabel('z');
        title('Induced charge on a grounded metal plane', 'FontWeight','bold');
        Q_ind_tot = sum(sum(Q_ind.*DIMENSIONSY/GRIDSIZEY).*DIMENSIONSZ/GRIDSIZEZ);
        set(gcf, 'Color', [1 1 1]);
       % hold on; % Add subsequent plots to the image
      %  surf(Q_ind);
      
      % Create a grid of x,y, and z values
      
      
    % plotting
        %plot the potential
%         phi_tot(phi_tot>0.5) = 0.5;
%         phi_tot(phi_tot<0) = NaN;
%         imagesc([-DIMENSIONSX/16 DIMENSIONSX/2],[-DIMENSIONSY/2 DIMENSIONSY/2], phi_tot(:, floor(7/16*GRIDSIZEX):GRIDSIZEX));  %plot the potential :, floor(GRIDSIZEX/2):GRIDSIZEX
%         set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','normal','ydir','normal','dataAspectRatio',[1 1 1]);     colorbar;
%         xlabel('x');
%         ylabel('y');
%         title('Induction on a grounded metal plane', 'FontWeight','bold');
%         hold on; % Add subsequent plots to    the image
% 
%         % change the E field arrays to plot less arrows in the field plot
%             arrowsarraydx = floor(size(xx,1)/NFIELDPLOTARROWS);
%             arrowsarraydy = floor(size(yy,1)/NFIELDPLOTARROWS);
%             Ex_quiver = Ex_tot(1:arrowsarraydx:end, 1:arrowsarraydy:end);
%             Ey_quiver = Ey_tot(1:arrowsarraydx:end, 1:arrowsarraydy:end);
%             xquiver = xx(1:arrowsarraydx:end, 1:arrowsarraydy:end);
%             yquiver = yy(1:arrowsarraydx:end, 1:arrowsarraydy:end);
%             Eabs_totquiver = Eabs_tot(1:arrowsarraydx:end, 1:arrowsarraydy:end);
%             Ex_quiver(Ex_quiver>1) = NaN;
%             Ex_quiver(Ex_quiver<-1) = NaN;
%             Ex_quiver(xquiver<0) = NaN;
%             Ey_quiver(Ey_quiver>1) = NaN;
%             Ey_quiver(Ey_quiver<-1) = NaN;
%             
%         % plot the E-Field
%         quiver(xquiver, yquiver, Ex_quiver./Eabs_totquiver, Ey_quiver./Eabs_totquiver,SIZEPLOTARROWS, 'Color','r');
%         fill([-0.1*DIMENSIONSX -0.1*DIMENSIONSX 0 0], [-DIMENSIONSY DIMENSIONSY DIMENSIONSY -DIMENSIONSY], 'white'); %rectangle to show the electrode 
%         fill([X0+0.005*DIMENSIONSX X0+0.005*DIMENSIONSX X0-0.002*DIMENSIONSX X0-0.002*DIMENSIONSX], [-0.001*DIMENSIONSY 0.001*DIMENSIONSY 0.001*DIMENSIONSY -0.001*DIMENSIONSY], 'black'); %line to show the electron 
%         fill([-0.1*DIMENSIONSX -0.1*DIMENSIONSX 0 0], [-0.001*DIMENSIONSY-3 0.001*DIMENSIONSY-3 0.001*DIMENSIONSY-3 -0.001*DIMENSIONSY-3], 'black'); %line to show the electrodes
%         fill([-0.1*DIMENSIONSX -0.1*DIMENSIONSX 0 0], [-0.001*DIMENSIONSY+3 0.001*DIMENSIONSY+3 0.001*DIMENSIONSY+3 -0.001*DIMENSIONSY+3], 'black'); %line to show the electrodes
%         hold off; % Add subsequent plots to the image
%         set(gcf, 'Color', [1 1 1]);
%         
% %         plot the E-Field at x = 0 in y direction
%         figure;
%         plot(floor(DIMENSIONSY/5):DIMENSIONSY/GRIDSIZEY:floor(DIMENSIONSY*4/5), Ex_tot(floor(size(Ex_tot,2)/5):floor(size(Ex_tot,2)*4/5),floor(size(Ex_tot,2)/2)+1), 'LineWidth',2, 'Color','red');      %plot Ex at x = 0
%         title('E-Field on the metal plane width charge at ', 'FontWeight','bold');
%         title_str = sprintf('E-Field on the metal plane width charge at x=%f', X0);
%         title(title_str, 'FontWeight','bold');
%         xlabel('y position [a.u.]');
%         ylabel('Electric field in x');
%         ylim([-1e-2 0]);
%         set(gcf, 'Color', [1 1 1]);   
 end %InducedCurrentMetallPlate  



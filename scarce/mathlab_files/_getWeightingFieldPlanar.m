function [E_x E_y] = getWeightingFieldPlanar(x,y,D,S)
    %from Nuclear Instruments and Methods in Physics Research A 535 (2004)
    %554–557, with correction from wbar = pi*w/2/D to wbar = pi*w/D
    
    %with x [um] is the position in the sensor, y [um] the offset from the
    %middle of the electrode, D [um] the sensor thickness and S [um] the
    %eletrode width. The field is calculated from the drivation of the potential
    %in x and y.
    
    y = D-y;    %electrode at D not at 0
    xbar = pi.*x./D;
    ybar = pi.*(y-D)./D;
    wbar = pi*S./D;
%     E_x = -1/pi.*( pi.*sin(ybar)./(2.*D.*(cosh(wbar./2+xbar)+cos(ybar))) -  pi.*sin(ybar)./(2.*D.*(cos(ybar)+cosh(wbar./2-xbar))) );    %not easy to find a more simple form
%     E_y = -1/pi.*( pi.*sinh(wbar./2+xbar)./(2.*D.*(cosh(wbar./2+xbar)+cos(1-pi.*y./D))) + pi.*sinh(wbar./2-xbar)./(2.*D.*(cos(1-pi.*y./D)+cosh(wbar./2-xbar))) ); %not easy to find a more simple form

    E_x = -sin(ybar)./(2.*D).*( 1./(cosh(xbar-wbar./2)+cos(ybar)) - 1./(cosh(wbar./2+xbar)+cos(ybar)) );
    E_y = -1./(2.*D).*( sinh(wbar./2-xbar)./(cosh(wbar./2-xbar)+cos(ybar)) + sinh(wbar./2+xbar)./(cosh(wbar./2+xbar)+cos(ybar)) );
    
end %getWeightingFieldPlanar
    
   



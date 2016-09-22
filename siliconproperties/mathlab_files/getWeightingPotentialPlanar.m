function Phi_w = getWeightingPotentialPlanar(x,y,D,S)
    %from Nuclear Instruments and Methods in Physics Research A 535 (2004)
    %554–557, with correction from wbar = pi*w/2/D to wbar = pi*w/D
    
    %with x [um] is the offset from the middle of the electrode, y [um] 
    %the position in the sensor, D [um] the sensor thickness and S [um] the
    %eletrode width.
    
    %wheighting potential for one pixel
    y = D-y;    %electrode at D not at 0
    xbar = pi.*x./D;
    ybar = pi.*(y-D)./D;
    STRIPWIDTHbar = pi*S./D;   
    Phi_w = -1./pi .* ( atan( tan(ybar./2).*tanh((xbar+STRIPWIDTHbar./2)./2) ) - atan( tan(ybar./2).*tanh((xbar-STRIPWIDTHbar./2)./2) ) );
end %getWeightingPotentialPlanar
    
   



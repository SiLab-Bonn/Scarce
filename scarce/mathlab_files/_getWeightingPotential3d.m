function Phi_w = getWeightingPotential3d(x,y,D,R)
    %from Ramo, Whinnery and van Duzer: Field and waves in communication
    %electronics p. 188
    
    %wheighting potential for two cylinders with radius R and distance D
    D = D/2;    %D is the total distance between the columns
    a = sqrt(D.*D-R.*R);
    Phi_w = 1./(4.*acosh(D./R)).*log( ((x-a).^2+y.^2)./((x+a).^2+y.^2) ) + 0.5;
    Phi_w(Phi_w<0) = NaN;
    Phi_w(Phi_w>1) = NaN;
    Phi_w(sqrt( (x+D).*(x+D) + y.*y )<R) = NaN;
    Phi_w(sqrt( (x-D).*(x-D) + y.*y )<R) = NaN;
end %getWeightingPotential3d
    
   



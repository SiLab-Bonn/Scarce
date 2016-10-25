function [E_x E_y] = getWeightingField3d(x,y,D,R)
    %from the analytical derivation of the getWeightingPotential function
    
    D = D/2;
    a = sqrt(D.*D-R.*R);
    
    E_x = a./(acosh(D./R)).*(a.^2-x.^2+y.^2)./( ((a-x).^2+y.^2) .* ((a+x).^2+y.^2) );
    E_y = -2.*a./(acosh(D./R)).*(x.*y)./( ((a-x).^2+y.^2) .* ((a+x).^2+y.^2) );
    
    E_x(sqrt( (x+D).*(x+D) + y.*y )<R) = NaN;
    E_x(sqrt( (x-D).*(x-D) + y.*y )<R) = NaN;
    E_y(sqrt( (x+D).*(x+D) + y.*y )<R) = NaN;
    E_y(sqrt( (x-D).*(x-D) + y.*y )<R) = NaN;
end %getWeightingField3d
    
   



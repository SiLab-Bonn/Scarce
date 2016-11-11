function [Q_tot, I_tot] = animateICUMP(startpoint, iterations)
    Q_tot = zeros(1, iterations);
    I_tot = zeros(1, iterations);
    Q_tot(1) = InducedCurrentMetallPlate(startpoint);
    I_tot(1) = Q_tot(1) - 0;
    f = getframe(gcf);
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,floor(iterations/10)) = 0;
    for j=1:iterations
        j
        Q_tot(1+j) = InducedCurrentMetallPlate(startpoint-j*startpoint/iterations);
        I_tot(1+j) = Q_tot(1+j) - Q_tot(j);
        f = getframe(gcf);
        im(:,:,1,j+1) = rgb2ind(f.cdata,map,'nodither');
        M(j) = getframe;
    end
    imwrite(im, map,'ICUMP.gif','DelayTime',0.010,'LoopCount',inf); %g443800
end %animateICUMP
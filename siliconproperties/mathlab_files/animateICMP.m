function animateICMP(startpoint, iterations)
    InducedChargeMetallPlate(startpoint);
    f = getframe(gcf);
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,floor(iterations/10)) = 0;
    for j=1:iterations
        InducedChargeMetallPlate(startpoint-j*startpoint/iterations);
        f = getframe(gcf);
        im(:,:,1,j+1) = rgb2ind(f.cdata,map,'nodither');
        M(j) = getframe;
    end
    imwrite(im, map,'ICMP.gif','DelayTime',0.010,'LoopCount',inf); %g443800
end %animateICMP
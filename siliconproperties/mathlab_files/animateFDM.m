function animateFDM(n)
    for j=1:100
        FDM(200,50,floor(1*j))
        M(j) = getframe;
    end
    %for j=110:20
    %    FDM(200,50,floor(10*j))
    %    M(j) = getframe;
    %end
    movie(M, 1, 4);
    movie2avi(M, 'test.avi', 'compression', 'none');
end %animateFDM
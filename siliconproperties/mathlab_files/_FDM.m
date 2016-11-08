%script to solve the laplace equation d/dx^2 + d/dy^2 = 0 with finite differences
%used to plot a potential distribution
function FDM(NUMBEROFPOINTS_X, NUMBEROFPOINTS_Y, ITERATIONS, CREATEPLOT)
    POTENTIALLEFT = 5;
    POTENTIALRIGHT = 4;
    POTENTIALTOP = 2;
    POTENTIALBOTTOM = 1;
    potential = zeros(NUMBEROFPOINTS_Y,NUMBEROFPOINTS_X);
    potential(:,1) = POTENTIALLEFT;
    potential(:,end) = POTENTIALRIGHT;
    potential(end,:) = POTENTIALTOP;
    potential(1,floor(NUMBEROFPOINTS_X./2-NUMBEROFPOINTS_Y./4):floor(NUMBEROFPOINTS_X./2+NUMBEROFPOINTS_Y./4)) = POTENTIALBOTTOM;
    if(CREATEPLOT)
        imagesc([-2 2], [0 1], potential);
        set(gcf, 'Position', [100 100 100+NUMBEROFPOINTS_X*4 100+NUMBEROFPOINTS_Y*4], 'Color', [1 1 1]);
        set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','normal','ydir','normal','dataAspectRatio',[1 1 1]);
        title('Iterative solving of the laplace equation with finite difference method, iteration 0', 'FontWeight','bold');
        %set(gca,'nextplot','replacechildren','visible','off')
        colorbar;
        f = getframe(gcf);
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,floor(ITERATIONS/10)) = 0;
    end;
    for(iter = 1:ITERATIONS)
        for(i_x = 2:NUMBEROFPOINTS_X-1)
            for(i_y = 2:NUMBEROFPOINTS_Y-1)
                potential(i_y,i_x) = (potential(i_y+1, i_x) + potential(i_y-1, i_x) + potential(i_y, i_x+1) + potential(i_y, i_x-1))./4;
            end
        end;
        if(CREATEPLOT && mod(iter,10) == 0)
           imagesc([-2 2], [0 1], potential);
           set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','normal','ydir','normal','dataAspectRatio',[1 1 1]);
           title_str = sprintf('Iterative solving of the laplace equation with finite difference method, iteration %d', iter);
           title(title_str, 'FontWeight','bold');
           colorbar;
           f = getframe(gcf);
           im(:,:,1,floor(iter./10)+1) = rgb2ind(f.cdata,map,'nodither');
           iter
        end;
    end
    %[Efield_x Efield_y] = gradient(-potential);
    %h = imagesc(potential);
    %set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','normal','ydir','normal','dataAspectRatio',[1 1 1]);
    %colorbar;
    
    %movie(M, 1, 4);
    if(CREATEPLOT)
        imwrite(im, map,'FDM.gif','DelayTime',0.010,'LoopCount',inf) %g443800
    end;
    %movie2avi(M, 'test.avi', 'compression', 'none');
    %eval(['print -djpeg99 '  num2str(ITERATIONS)]);
    %ITERATIONS
    %imwrite(potential, 'test.jpg', 'jpg');
    %hold on;
    %quiver(Efield_x,Efield_y, 2, 'Color','r');
    %hold off;
end %FDM
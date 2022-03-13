directionStr='';

numConcavitiesTested=30;
colormp=copper(numConcavitiesTested)

for i=1:numConcavitiesTested
        %currConcavityFrac=(i-1-numConcavitiesTested/2+0.5)/(numConcavitiesTested/2); %from -1 to 1
         %currConcavityFrac=-(i-1-numConcavitiesTested/2+0.5)/(numConcavitiesTested/2);
         currConcavityFrac=((i-1-numConcavitiesTested/2+0.5)/(numConcavitiesTested/2));

    %[curveX,curveY] = drawBezier([0 0 currConcavityFrac 1],[1 currConcavityFrac 0 0],100);
    if(currConcavityFrac<0)
        %[curveX,curveY] = drawBezier([0 -currConcavityFrac 1 1],[1 1 -currConcavityFrac 0],100);
           [curveX,curveY] = drawBezier([0 0 1+currConcavityFrac  1],[0 -currConcavityFrac 1 1],100);

    else
        %[curveX,curveY] = drawBezier([0 0 1-currConcavityFrac 1],[1 1-currConcavityFrac 0 0],100);
           [curveX,curveY] = drawBezier([0 currConcavityFrac 1 1],[0 0 1-currConcavityFrac 1],100);

    end
    if(strcmp(directionStr,'proxToDist'))
        plot((1-curveX),curveY,'Color',colormp(i,:),'LineWidth',3)
    else
       plot(curveX,curveY,'Color',colormp(i,:),'LineWidth',3)
    end
    hold on
    %ylabel('Relative spike time (msec)')
    %xlabel('Distance along dendrite from distal end (um)')

end
title('Varying convexity with Bezier curves')
colormap(gca,copper)
    cb=colorbar
    ylabel(cb,'Convexity')
    caxis([-1 1])
function [cmiDataVal,cmiBootstrapHistCounts,cmiBootstrapHistBinCenters] = getConditionalMutualInfo(x,y,z,xStr,yStr,zStr)
%given 3 r.v. observations column vectors returns mutual info of first two
%conditional on third
showPlots=0;

numShuffles=500;

varNames={'Time','Speed','Distance','Phase'};

for vi=1:length(varNames)
    currVarName=varNames{vi};
    
    if(contains(xStr,currVarName))
        xSaveStr=currVarName;
    end

    if(contains(zStr,currVarName))
        zSaveStr=currVarName;
    end


    if(contains(yStr,currVarName))
        ySaveStr=currVarName;
    end
end

%{
if(~exist('xStr'))
    xStr='Time (field frac)';
end


if(contains(xStr,'Time'))
    xSaveStr='time';
else
     xSaveStr='speed';
end
    yStr='Phase (cycle frac)';
    zStr='Distance (field frac)';
%}

cmiOfShuffledData=NaN(numShuffles,1);
%first loop is real data, all others are shuffles
for shi=1:(numShuffles+1)
    if(mod(shi,100)==0)
        disp(shi/numShuffles)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %assumes all in range [0, 1] (t,x,p)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %x=x/max(x(:));
    %y=y/max(y(:));
    %z=z/max(z(:));
    

    if(shi==1)
        x(x>1)=NaN;
        y(y>1)=NaN;
        z(z>1)=NaN;
        
        anyNaNrow=(isnan(x) | isnan(y) | isnan(z));
        x(anyNaNrow)=[];
        y(anyNaNrow)=[];
        z(anyNaNrow)=[];
        
        showPlots=1;
    else
        showPlots=0;
    end
    

    if(shi>1)
         %x=x(randperm(length(x)));
        y=y(randperm(length(y)));
         z=z(randperm(length(z)));
    end

    commonBinWidth=0.025;
    commonNumBins=1/commonBinWidth; %change in histcn code...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %prepare prerequisite distributions (smoothed)
    %I(X,Y|Z) = sumZ(sumY(sumX(pXYZ*log2(pZ*pXYZ/(pXZ*pYZ)))))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pZn pZEdges pZMids zInBins] = histcn(z(:));
    pZn=pZn(1:(end-1));
    pZ=pZn/sum(pZn(:));
   
       %univSmoothParam=0.05;
      univSmoothParam=0.15;
      
       %univSmoothParam=0.05;
        %  univSmoothParam=0.02;

    pZSmoothed=smooth(pZ,univSmoothParam); %default 5 bin local smoothing


    [pXZn pXZEdges pXZMids xzInBins] = histcn([x(:) z(:)]);
    pXZn=pXZn(1:(end-1),1:(end-1));
    pXZ=pXZn/sum(pXZn(:));
    defaultSmoothParamXZ=0.0111;
    [pXZSmoothed,smoothParamXZ]=smoothn(pXZ,univSmoothParam); 
    [pXZSmoothed,smoothParamXZ]=smoothn(pXZSmoothed,univSmoothParam); 


    [pYZn pYZEdges pYZMids yzInBins] = histcn([y(:) z(:)]);
    pYZn=pYZn(1:(end-1),1:(end-1));
    pYZ=pYZn/sum(pYZn(:));
    defaultSmoothParamYZ=0.0308;
    [pYZSmoothed,smoothParamYZ]=smoothn(pYZ,univSmoothParam); 
    [pYZSmoothed,smoothParamYZ]=smoothn(pYZSmoothed,univSmoothParam); 


    [pXYZn pXYZEdges pXYZMids xyzInBins] = histcn([x(:) y(:) z(:)]);
    pXYZn=pXYZn(1:(end-1),1:(end-1),1:(end-1));
    pXYZ=pXYZn/sum(pXYZn(:));
    defaultSmoothParamXYZ=0.0527;
    [pXYZSmoothed,smoothParamXYZ]=smoothn(pXYZ,univSmoothParam); 
    [pXYZSmoothed,smoothParamXYZ]=smoothn(pXYZSmoothed,univSmoothParam); 
    [pXYZSmoothed,smoothParamXYZ]=smoothn(pXYZSmoothed,univSmoothParam); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot to check for smoothness (representative of underlying probability
    %distributions?)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(showPlots)
        fH=figure
        subplot(2,2,1)
        plot(pZMids{1},pZSmoothed,'k-','LineWidth',3)
        hold on
        plot(pZMids{1},pZSmoothed,'ko')
        xlabel(zStr)
        ylabel('probability')
        title('pZ')

        subplot(2,2,2)
        omarPcolor(pXZMids{1},pXZMids{2},pXZSmoothed',fH)
        colormap(gca,jet)
        cb1=colorbar
        ylabel(cb1,'probability')
        xlabel(xStr)
        ylabel(zStr)
        title('pXZ')
        daspect([1 1 1])

        subplot(2,2,3)
        omarPcolor(pYZMids{2},pYZMids{1},pYZSmoothed,fH)
        colormap(gca,jet)
        cb2=colorbar
        ylabel(cb2,'probability')
        xlabel(zStr)
        ylabel(yStr)
        title('pYZ')
        daspect([1 1 1])

        subplot(2,2,4)
        %scatter3(pXYZMids{1},pXYZMids{2},pXYZMids{3}, squeeze(pXYZ(1,1,:)), squeeze(pXYZ(1,:,:)), 'filled')
        %set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
        zSlices=linspace(0,max(z(:))-0.01,32);
        dispPxyz=pXYZ;
        threshDispP=prctile(pXYZ(:),90);
        dispPxyz(dispPxyz<threshDispP)=NaN;
        for zi=1:length(zSlices)
            s=slice(pXYZMids{1},pXYZMids{2},pXYZMids{3}, dispPxyz, 0, 0, zSlices(zi));
            if(zi==1)
                xlabel(xStr)
                ylabel(yStr)
                zlabel(zStr)
                colormap(gca,jet)
                cb3=colorbar
                ylabel(cb3,'probability')
                hold on
                    daspect([1 1 1])
                    title('pXYZ')
                    maxFig
                    setFigFontTo(18)
            end  
        set(s, 'EdgeColor', 'none')
            drawnow
            pause(0.2)
            alpha 0.3
        end
        
        saveas(gcf,sprintf('%sVs%sGiven%sDistributions.tif',xSaveStr,ySaveStr,zSaveStr))
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute conditional mutual information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    distrLength=length(pZ);

    totalDevFromCondIndep=0;
    for zi=1:distrLength
        for yi=1:distrLength
            for xi=1:distrLength
                %currProbWeight=pXYZ(xi,yi,zi);
                %currProbWeight=pXYZSmoothed(xi,yi,zi);

                %numer=currProbWeight*pZ(zi);
                numer=pXYZSmoothed(xi,yi,zi)*pZSmoothed(zi);
                %denom=pXZ(xi,zi)*pYZ(yi,zi);
                 denom=pXZSmoothed(xi,zi)*pYZSmoothed(yi,zi);

                if(numer<=0 || denom<=0)
                    continue; %convention is log2(0) is 0, add nothing, 
                     %also skips points where independence predicts zero, out of range
                end

                currDevFromCondIndep=pXYZSmoothed(xi,yi,zi)*log2(numer/denom);

                totalDevFromCondIndep=totalDevFromCondIndep+currDevFromCondIndep;
            end
        end

    end

    currCMIVal=totalDevFromCondIndep;
    if(shi==1)
        cmiDataVal=currCMIVal;
    elseif(shi>1)
        cmiOfShuffledData(shi-1)=currCMIVal;
    end

end

figure; 
hH=histogram(cmiOfShuffledData);
hold on
pH=plot([cmiDataVal cmiDataVal],ylim,'k--','LineWidth',5)

xlabel('CMI (deviation from conditional independence)')

ylabel('Count')
legend([hH pH],{'shuffled data cmi','original data cmi'})
title({sprintf('%s Vs %s, controlling for %s',xSaveStr,ySaveStr,zSaveStr),sprintf('Bootstrap CMI distribution, %d shuffles',numShuffles)})

maxFig
setFigFontTo(18)
saveas(gcf,sprintf('%sVs%sGiven%sInfoTheoryStats.tif',xSaveStr,ySaveStr,zSaveStr))
disp('')


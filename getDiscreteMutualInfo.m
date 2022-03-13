function [Ixy,smoothExcessProbability,Hx,Hy,rhoSpear,pSpear,dispXlims,dispYlims] = getDiscreteMutualInfo(x,y,nbinsX,nbinsY,xLabelName,yLabelName,showPlots)
    %minEstimableProb=5e-3;
     %minEstimableProb=1e-3;
     %minEstimableProb=5e-4;
     %setTightSubplots_SpaceTime
     minEstimableProb=0;
     plotDeviationFromInd=0;
     
     %showPlots=1;
     %showPlots=0;
     fitLineToMutualDep=1;
    
    [xInBinNums,xEdges]=discretize(x,nbinsX);
    [yInBinNums,yEdges]=discretize(y,nbinsY);
    
    dispXlims=[min(xEdges) max(xEdges)];
    dispYlims=[min(yEdges) max(yEdges)];
    
    xBinCenters=edgesToBins(xEdges);
    yBinCenters=edgesToBins(yEdges);
    
    if(showPlots)
        fH=figure;
        subplot(2,2,3)
        h=histogram2(x,y,xEdges,yEdges, 'DisplayStyle','tile');
        colormap(jet)
        xlim(dispXlims)
        ylim(dispYlims)
        jointPMF=h.Values;
    else
        [jointPMF,~,~]=histcounts2(x,y,xEdges,yEdges);
    end
    
    
    %jointPMF=imgaussfilt(jointPMF,1);
    
    jointPMF=jointPMF/sum(jointPMF(:));
    [xMarginal,~]=histcounts(x,xEdges);
    [yMarginal,~]=histcounts(y,yEdges);
    
    xMarginal=xMarginal/sum(xMarginal(:));
    yMarginal=yMarginal/sum(yMarginal(:));
    Hx=getEntropy(xMarginal);
    Hy=getEntropy(yMarginal);
    
    if(showPlots)
        %subplot(2,2,1)
        subplot(3,2,1)
        plot(xBinCenters,xMarginal,'k-')

        title(sprintf('H(X)=%.5f',Hx))
        xlim(dispXlims)
        ylim([0 Inf])
        xlabel(xLabelName)
        ylabel('Probability')

        subplot(3,2,4)
        plot(yMarginal,yBinCenters,'k-')
         ylabel(yLabelName)
         xlabel('Probability')
         title(sprintf('H(Y)=%.5f',Hy))
        xlim([0 Inf])
        ylim(dispYlims)
    end
    


    Ixy=0;
    %dataJointProb=NaN(size(y,1)*size(x,1),1);
    %multModelJointProb=NaN(size(y,1)*size(x,1),1);
    dataJointProb=zeros(nbinsX*nbinsY,1);
    multModelJointProb=zeros(nbinsX*nbinsY,1);
    
    count=0;
    %loop over value in each bin, not each data point
    %for yi=1:size(y,1)
    %    for xi=1:size(x,1)
    for currY=min(yInBinNums):max(yInBinNums)
            %currX=xInBinNums(xi);
            %currY=yInBinNums(yi);
         
        for currX=min(xInBinNums):max(xInBinNums)
            
            count=count+1;
          
            if(xMarginal(currX)*yMarginal(currY)>minEstimableProb && jointPMF(currX,currY)>minEstimableProb)
                dataJointProb(count)=jointPMF(currX,currY);
                multModelJointProb(count)=xMarginal(currX)*yMarginal(currY);
                Ixy=Ixy+dataJointProb(count)*log2(dataJointProb(count)/(multModelJointProb(count)));
            end
        end
    end %y bin loop
    dataJoint=reshape(dataJointProb,[nbinsX nbinsY]);
    multModelJoint=reshape(multModelJointProb,[nbinsX nbinsY]);
    multModelJoint(multModelJoint==0)=min(dataJoint(dataJoint>0)); %no zero in denominator

    %smoothExcessProbability=imgaussfilt(dataJoint'./multModelJoint',5,'Padding','circular')
    %smoothExcessProbability=imgaussfilt(dataJoint'./multModelJoint',3,'Padding','replicate');
    smoothExcessProbability=imgaussfilt(dataJoint'./multModelJoint',3,'Padding','replicate');
        
    if(fitLineToMutualDep)
        %figure
       overDenseImgLogical=smoothExcessProbability>1;
       overDenseImg=zeros(size(overDenseImgLogical));
       overDenseImg(overDenseImgLogical)=1;
        currColCOM=NaN(size(overDenseImg,2),1);

        for col=1:size(overDenseImg,2)
            %currCol=(overDenseImg(:,col)/nbinsY)*2*pi;
            currCol=(smoothExcessProbability(:,col)/nbinsY)*2*pi;
            if(max(currCol)==0)
                currColCOM(col)=NaN; %No preferred phase
            else
                currColCOM(col)=getCircCOM(currCol);
            end
        end
        
        phasesInRadians=(currColCOM)*2*pi/(max(currColCOM));
        %linearPos=(1:nbinsX)/nbinsX;
        
         %subplot(3,2,5)
        
        %omarPcolor(xBinCenters,yBinCenters,(overDenseImg),fH)
        %title('over concentrated joint density')
        %xlabel(xLabelName)
        %ylabel([yLabelName ' center of mass'])
        
        subplot(3,2,6)
        plot(xBinCenters,phasesInRadians)
	ylim([0 2*pi])
	xlim(dispXlims)
        [ rhoK,pK,sK,bK ] = kempter_lincirc(xBinCenters(:),phasesInRadians(:));
        
        [rhoSpearman, pSpearman]=corr([xBinCenters(:) phasesInRadians(:)],'Type','Spearman');
        
        rhoSpear=rhoSpearman(2,1);
        pSpear=pSpearman(2,1);
        
        xlabel(xLabelName)
        ylabel(['COM ' yLabelName ])
        
        title({sprintf('p_{corr}=%0.5f,slope=%.5f radians/field, rho=%.5f',pK,sK*360,rhoK),sprintf('CircCOM Spearman rank coeff: %.5f, p=%.5f',rhoSpear,pSpear)})
    end

    if(showPlots)
        %figure
        %subplot(2,2,3)
        subplot(3,2,3)
        
        %omarPcolor(xBinCenters,[yBinCenters yBinCenters+360] ,[smoothExcessProbability; smoothExcessProbability],fH)
        omarPcolor(xBinCenters,[yBinCenters] ,[smoothExcessProbability],fH)
        colormap(jet)
        cb=colorbar('north')
        %cb.Position=[0.55 0.75 0.4 0.05 ];
        cb.Position=[0.125 0.275 0.35 0.025 ];
        ylabel(cb,'Joint density relative to independence')
        caxis([0.8 1.3])
        %caxis([0.8 2])
            %subplot(2,2,3)
             subplot(3,2,3)
        title(sprintf('I(X;Y) = %0.5f',Ixy))
        setFigFontTo(18)
        xlabel(xLabelName)
        ylabel(yLabelName)
        %title('Joint dist. assuming independence')
        %set(gca,'xscale','log')
        %xlim([0 1000])
        %ylim([0 360])
        %set(gca,'yscale','log')
        
    end
    
    if(plotDeviationFromInd)
        figure;
        %subplot(1,2,1)
        plot(multModelJointProb,dataJointProb,'r.')
        hold on
        maxPoint=max(multModelJointProb(:));
        plot([0 maxPoint] ,[0 maxPoint],'k-','LineWidth',3)
        ylabel('Joint probability in data')
        xlabel('Joint probability in model of independence')
        title('Joint probability deviation from independence')
        xlim([-Inf Inf])
        ylim([-Inf Inf])

        daspect([1 1 1])

        %subplot(1,2,2)
        %jointToModelRatio=dataJointProb./multModelJointProb;
        %histogram(jointToModelRatio)
        setFigFontTo(18)
        %disp('here')
    end
    

function [] = plotCrossSections(xCenters,yCenters,zValues,xBinWidth,yBinWidth,xLabel,yLabel,zLabel)
%given a surface, plots cross sections for each bin in x and y dimensions

    holdXEdges=xCenters(1):xBinWidth:xCenters(end);
    holdXBinCenters=edgesToBins(holdXEdges);
    
    holdYEdges=yCenters(1):yBinWidth:yCenters(end);
    holdYBinCenters=edgesToBins(holdYEdges);
    
    crossSectionRes=50;
       yQVals=linspace(yCenters(1),yCenters(end),crossSectionRes);
       xQVals=linspace(xCenters(1),xCenters(end),crossSectionRes);
       
       
       xColors=jet(length(holdXBinCenters));
       yColors=jet(length(holdYBinCenters)); 
    
    figure
    subplot(1,2,1)
    %plot holding x
    zCurveAtGivenX=[];
    for i=1:length(holdXBinCenters)
        currHoldXVal=holdXBinCenters(i);
     
        zCurveAtGivenX{i}=interp2(xCenters,yCenters,zValues,currHoldXVal,yQVals);
        plot(yQVals,zCurveAtGivenX{i},'-','LineWidth',4,'Color',xColors(i,:))
        hold on
    end
    
    cb=colorbar
    colormap(jet)
    xlabel(yLabel)
    ylabel(zLabel)
    ylabel(cb,xLabel)
    caxis([min(xCenters) max(xCenters)])
    
    
    subplot(1,2,2)
    %plot holding x
    zCurveAtGivenY=[];
    for i=1:length(holdYBinCenters)
        currHoldYVal=holdYBinCenters(i);
    
        zCurveAtGivenY{i}=interp2(xCenters,yCenters,zValues,xQVals,currHoldYVal);
        plot(xQVals,zCurveAtGivenY{i},'-','LineWidth',4,'Color',yColors(i,:))
        hold on
    end
     xlabel(xLabel)
    ylabel(zLabel)
    
     cb=colorbar;
    colormap(jet)
    ylabel(cb,yLabel)
    
    caxis([min(yCenters) max(yCenters)])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot conditional change densities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phaseChangeMin=-100; %deg
    phaseChangeMax=100; %deg
    figure
    
    allGivenXChangeCoords=[];
    allGivenYChangeCoords=[];
    
    allGivenX_ZYslopes=[];
    allGivenY_ZXslopes=[];
    
    
    for i=1:length(holdXBinCenters)
        currZYgivenX = zCurveAtGivenX{i};
        currZXgivenY = zCurveAtGivenY{i};
        
        validIdxesGivenX=~isnan(currZYgivenX);
        validIdxesGivenY=~isnan(currZXgivenY);
        
        currZYgivenX=currZYgivenX(validIdxesGivenX);
        currZXgivenY=currZXgivenY(validIdxesGivenY);

        
        currZYgivenXzeroed=currZYgivenX-nanmean(currZYgivenX(:));
        currZXgivenYzeroed=currZXgivenY-nanmean(currZXgivenY(:));
        
        currXqVals=xQVals(validIdxesGivenY);
        currYqVals=yQVals(validIdxesGivenX);
        
        
        currXqValsZeroed=currXqVals-nanmean(currXqVals);
        currYqValsZeroed=currYqVals-nanmean(currYqVals);
        
        allGivenXChangeCoords=[allGivenXChangeCoords; [currYqValsZeroed(:), currZYgivenXzeroed(:)]];

        allGivenYChangeCoords=[allGivenYChangeCoords; [currXqValsZeroed(:), currZXgivenYzeroed(:)]];
        
        [mZY,b,R,p]=getLinearFit(currYqValsZeroed,currZYgivenXzeroed);
        
        [mZX,b,R,p]=getLinearFit(currXqValsZeroed,currZXgivenYzeroed);
        
        allGivenX_ZYslopes(i)=mZY;
        allGivenY_ZXslopes(i)=mZX;
        
            subplot(1,2,2)
            plot(currYqValsZeroed,currZYgivenXzeroed,'k-','LineWidth',1)
            hold on
            ylim([phaseChangeMin phaseChangeMax])
        
            axis square
            
            subplot(1,2,1)
            
            plot(currXqValsZeroed,currZXgivenYzeroed,'k-','LineWidth',1)
            hold on
            ylim([phaseChangeMin phaseChangeMax])
            
             axis square
        
    end
    
    subplot(1,2,2)
    vline(0,'k',1)
    hline(0,'k',1)
    box off
     xlabel('\DeltaTime')
    ylabel('\DeltaPhase')
    subplot(1,2,1)
    vline(0,'k',1)
    hline(0,'k',1)
    box off
    xlabel('\DeltaDistance')
    ylabel('\DeltaPhase')
    %{
    fChangeH=figure;
    numBins=100;
    
    xChangeEdges=linspace(-0.4,0.4,numBins+1);
    yChangeEdges=linspace(-0.4,0.4,numBins+1);
    phaseChangeEdges=linspace(phaseChangeMin,phaseChangeMax,numBins+1);
    
    subplot(1,2,1)
 getJointDistr(allGivenXChangeCoords(:,1),allGivenXChangeCoords(:,2),yChangeEdges,phaseChangeEdges,fChangeH);
     ylim([phaseChangeMin phaseChangeMax])
     
    subplot(1,2,2)
     getJointDistr(allGivenYChangeCoords(:,1),allGivenYChangeCoords(:,2),xChangeEdges,phaseChangeEdges,fChangeH);

 ylim([phaseChangeMin phaseChangeMax])
    %}
    
    disp('')
    
    figure
    zyMean=nanmean(allGivenX_ZYslopes);
    zySEM=getSEMacrossRows(allGivenX_ZYslopes(:));
    
    zxMean=nanmean(allGivenY_ZXslopes);
    zxSEM=getSEMacrossRows(allGivenY_ZXslopes(:));
    
    bar(abs([zxMean zyMean ]),'FaceColor','none','LineWidth',3)
    hold on
    errorbar([1 2],abs([zxMean zyMean ]),[zxSEM zySEM ],'k.')
    
    ylim([0 270])
    box off
    
    [p,h,stats] = ranksum(allGivenY_ZXslopes,allGivenX_ZYslopes)
    
    disp('')
    
    
    

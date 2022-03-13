function [devFromLinearVsSpace] = getLinearDeviationOfCylData(phasesDeg,distNorm,distBinWidth,titleStr)
%given a list of phases (0 to 360) and linear distances (between 0 and 1)
%returns the avg. deviation from linear relationship between phase and distance
%on 3D cylinder for each distance bin from 0 to 1
    close all
    
    clearvars -except phasesDeg distNorm distBinWidth titleStr
    findBestLinearCorrShift=1;
    if(findBestLinearCorrShift)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %circularly shift phases until error from model linear precession is
        %minimal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        degShifts=1:360;
        rhos=NaN(size(degShifts));
        totalVars=NaN(size(degShifts));
        for si=1:length(degShifts)
            currShift=degShifts(si);
            shiftedPhases=mod(phasesDeg+currShift,360);
            nonNaNids=~isnan(shiftedPhases);
            totalVars(si)=var(shiftedPhases(nonNaNids));
            currRho=corr(distNorm(nonNaNids)',shiftedPhases(nonNaNids));
            rhos(si)=currRho;
        end
        [~,maxRhoID]=max(abs(rhos));
        %[~,minVarID]=min(totalVars);
        bestPhaseShift=degShifts(maxRhoID)
        %bestPhaseShift=degShifts(minVarID)
        %{
        figure; 
        yyaxis left
        plot(totalVars)
        xlabel('degree shift')
        ylabel('Total variance')
        
        yyaxis right
        plot(rhos)
        xlabel('degree shift')
        ylabel('Linear correlation')
        %}

        phasesDeg=mod(phasesDeg+bestPhaseShift,360);
    end
    
    figure;

    plot(distNorm,phasesDeg,'k.')
    hold on
    linearPhaseModel=(1-(distNorm(:)))*360;
    plot(distNorm,linearPhaseModel,'b-','LineWidth',3)
    legend('data','linear null model')
    ylim([0 360])
    xlim([0 1])
    xlabel('Distance in field (normalized)')
    ylabel('Theta Phase')
    setFigFontTo(18)
    
    saveas(gcf,sprintf('%s.tif',titleStr))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %convert input into (x,y,z) triplets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phasesRadians=phasesDeg*pi/180;

    xCoords=cos(phasesRadians(:));
    yCoords=sin(phasesRadians(:));
    zCoords=distNorm(:);

    numPts=length(zCoords);
    
    pts=[xCoords,yCoords,zCoords];

    %error checking
    assert(numPts==length(xCoords) && numPts==length(yCoords),'number of phases and number of distances must be same')

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %create linear reference curve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numBins=1/distBinWidth;
    angularBinSize=2*pi/numBins;
    
    
    angularBinCenters=(angularBinSize/2):angularBinSize:(2*pi-angularBinSize/2);
    binCenters=(distBinWidth/2):distBinWidth:(1-distBinWidth/2);
    binEdges=linspace(0,1,numBins+1);
    
    linearRefAngles=(1-binCenters)*2*pi; %linear phase precession points
    linearRefAngleEdges=(1-binEdges)*2*pi;
    
    linearRefPts=[cos(linearRefAngles(:)),sin(linearRefAngles(:)),binCenters(:)];
    linearRefPtEdges=[cos(linearRefAngleEdges(:)),sin(linearRefAngleEdges(:)),binEdges(:)];
    
    linearRefTangents=NaN(numBins,3);
    linearRefAxials=NaN(numBins,3);
    linearRefPerps=NaN(numBins,3);
    
    linearRefTangentEdges=NaN(numBins+1,3);
    linearRefAxialEdges=NaN(numBins+1,3);
    linearRefPerpEdges=NaN(numBins+1,3);
    
    %figure
    for ri=1:(numBins) %bin center ri corresponds to lower edge for edge vars
        centerPt=[0 0 binCenters(ri)];
        centerPtLowerEdge=[0 0 (binCenters(ri) - distBinWidth/2)];
        centerPtUpperEdge=[0 0 (binCenters(ri) + distBinWidth/2)];
        
        if(ri+1>numBins)
            linearRefTangents(ri,:)=[0 0 0];
            linearRefAxials(ri,:)=[0 0 0];
        else
            linearRefTangents(ri,:)=linearRefPts(ri+1,:)-linearRefPts(ri,:);
            linearRefTangentEdges(ri,:)=linearRefPtEdges(ri+1,:)-linearRefPtEdges(ri,:);
            
            
            linearRefTangents(ri,:)=normVec(linearRefTangents(ri,:));
            linearRefTangentEdges(ri,:)=normVec(linearRefTangentEdges(ri,:));
            
            %towards center of circle is perpendicular (axial)
            linearRefAxials(ri,:)=(centerPt)-linearRefPts(ri,:); 
            linearRefAxialEdges(ri,:)=(centerPtLowerEdge)-linearRefPtEdges(ri,:); 
            
            %linearRefAxials(ri,3)=0; %should be flat in z

            linearRefAxials(ri,:)=normVec(linearRefAxials(ri,:));
            linearRefAxialEdges(ri,:)=normVec(linearRefAxialEdges(ri,:));
        end
        %linearRefTangents(ri,:)=linearRefPts(ri+1,:);
        
        linearRefPerps(ri,:)=cross(linearRefTangents(ri,:),linearRefAxials(ri,:));
        linearRefPerpEdges(ri,:)=cross(linearRefTangentEdges(ri,:),linearRefAxialEdges(ri,:));
       
          %plot3Dvectors([linearRefTangents(ri,:);linearRefPerps(ri,:)]/10,linearRefPts(ri,:),'kr')
          %plot3DvectorFromTo(linearRefPts(ri,:),linearRefPts(ri+1,:),'k')
          %plot3DvectorFromTo(linearRefPts(ri,:),centerPt+[0 0 0.5],'b')
          
          %plot3DvectorFromTo(linearRefTangents(ri,:))
          %plot3DvectorFromTo(linearRefPerps(ri,:))
    end
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop through bins
    %to store deviation from linear reference
    %for each distance bin
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    meanLinearDevPerBin=NaN(numBins,1);
    fCrossSections=figure;
    maxFig
    subplot(2,1,1)
    
    referenceColors=jet(numBins);
    for i=1:numPts
        currPtPlot=pts(i,:);
        plot3(currPtPlot(1),currPtPlot(2),currPtPlot(3),'.')
        hold on
    end
    
    for bi=1:numBins %bin center bi corresponds to lower edge for edge vars
        subplot(2,1,1)
        currBinCenter=binCenters(bi); 
        currBinRefPt=linearRefPts(bi,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find bounds perpedicular to control spiral
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %boundVectorNormalLower=-linearRefPerpEdges(bi,:);
        %boundVectorNormalUpper=-linearRefPerpEdges(bi+1,:);
        boundVectorNormalLower=linearRefTangentEdges(bi,:);
        boundPointLower=linearRefPtEdges(bi,:);
        %boundVectorNormalUpper=-linearRefTangentEdges(bi+1,:);
        boundVectorNormalUpper=-linearRefTangentEdges(bi,:);
        boundPointUpper=linearRefPtEdges(bi+1,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %each bound plane forms a half-space, must be on correct side of
        %each to be counted in this bin
        %third plane bound is diameter perpendicular to ref pt angle, any z
        %normal is point from current z to current bin center reference pt
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        boundVectorNormalAxial=currBinRefPt-[0 0 currBinCenter];
        boundPointAxial=[0 0 currBinCenter];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %fourth,fifth plane bound is parallel to tangent, starting at
        %z+-0.8 from current bin center reference pt (assumes reasonable range of
        %slopes)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zMaxSearch=0.6;
        boundVectorNormalDistanceUpper=linearRefPerpEdges(bi,:);
        boundPointDistanceUpper=currBinRefPt+[0 0 zMaxSearch];
        
        boundVectorNormalDistanceLower=-linearRefPerpEdges(bi,:);
        boundPointDistanceLower=currBinRefPt-[0 0 zMaxSearch];
        
        
        [isOnNormalSideLower]=pointIsOnNormalSide(boundPointLower,boundVectorNormalLower,pts);
        [isOnNormalSideUpper]=pointIsOnNormalSide(boundPointUpper,boundVectorNormalUpper,pts);
        [isOnNormalSideCenter]=pointIsOnNormalSide(boundPointAxial,boundVectorNormalAxial,pts);
        [isOnNormalSideDistanceLower]=pointIsOnNormalSide(boundPointDistanceLower,boundVectorNormalDistanceLower,pts);
        [isOnNormalSideDistanceUpper]=pointIsOnNormalSide(boundPointDistanceUpper,boundVectorNormalDistanceUpper,pts);
        
        isInCurrBin=isOnNormalSideLower & isOnNormalSideUpper & isOnNormalSideCenter & isOnNormalSideDistanceLower ...
                & isOnNormalSideDistanceUpper;
        
        
        currRefZ=currBinRefPt(3);
        
        currBinDataPts=pts(isInCurrBin,:);
        
        
        numPtsCurrBin=size(currBinDataPts,1);
        hold on
        
        scatter3(currBinDataPts(:,1),currBinDataPts(:,2),currBinDataPts(:,3),50*ones(numPtsCurrBin,1),referenceColors(bi,:))
        scatter3(currBinRefPt(:,1),currBinRefPt(:,2),currBinRefPt(:,3),100,'k')
        
        p1H=drawPlaneGivenNormal(boundPointLower,boundVectorNormalLower);
        p2H=drawPlaneGivenNormal(boundPointUpper,boundVectorNormalUpper);
        axialLims=[-1 1];
        distPlaneLimsX=[currBinRefPt(1)-0.5 currBinRefPt(1)+0.5];
        distPlaneLimsY=[currBinRefPt(2)-0.5 currBinRefPt(2)+0.5];
        p3H=drawPlaneGivenNormal(boundPointAxial,boundVectorNormalAxial,axialLims,axialLims);
        %p3H.Visible='off';
        %p4H=drawPlaneGivenNormal(boundPointDistanceUpper,boundVectorNormalDistanceUpper,distPlaneLimsX,distPlaneLimsY,1);
        %p5H=drawPlaneGivenNormal(boundPointDistanceLower,boundVectorNormalDistanceLower,distPlaneLimsX,distPlaneLimsY,1);
        
        p4H=drawPlaneGivenNormal(boundPointDistanceUpper,boundVectorNormalDistanceUpper);
        p5H=drawPlaneGivenNormal(boundPointDistanceLower,boundVectorNormalDistanceLower);
        daspect([1 1 1])
        xlim([-1 1])
        ylim([-1 1])
        zlim([0 1])
   
        
        %plot3DvectorFromTo([0 0 0], linearRefPerps(biFixed,:),'b')
        %plot3DvectorFromTo([0 0 0], linearRefTangents(biFixed,:),'k')
        %plot3DvectorFromTo([0 0 0], majorAxisVec,'b')
        %plot3DvectorFromTo([0 0 0], minorAxisVec,'k')
        %plot3DvectorFromTo(currBinRefPt, currBinRefPt-linearRefPerps(bi,:),'r')
        %plot3DvectorFromTo(currBinRefPt, currBinRefPt+linearRefAxials(bi,:),'r')
        
        %view(80,16)
        view(38.4,22.4)
        drawnow
        %pause(0.01)
        
        currBinPtDeviations=NaN(numPtsCurrBin,1);
        
        
        for pti=1:numPtsCurrBin
            currInBinPt=currBinDataPts(pti,:)';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %if dist less than ref, deviation is negative (sublinear vs. superlinear)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                currBinPtDeviations(pti)=pdist([currBinRefPt(:)'; currInBinPt(:)']);
                
            if(currInBinPt(3)<currRefZ)
                currBinPtDeviations(pti)=-currBinPtDeviations(pti);
            end
        end
        
        meanLinearDevPerBin(bi)=nanmean(currBinPtDeviations);
         %meanLinearDevPerBin(bi)=nanmedian(currBinPtDeviations);
        
        zlabel('Distance in field')
        xlabel('cos(Theta phase)')
        ylabel('sin(Theta phase)')
        
        title(removeUnderscores(titleStr))
        
        subplot(2,2,3)
     
        
        %plot(binCenters,meanLinearDevPerBin,'k-','LineWidth',3)
        cumDev=cumsum(meanLinearDevPerBin,'omitnan');
        pD=plot(binCenters(1:bi),cumDev(1:bi),'k-','LineWidth',3);
        hold on
         plot([0 1], [0 0],'k--','LineWidth',3)
        ylabel('Total Deviation from null linear model')
        xlabel('Distance along curve')
        xlim([0 1])
        %ylim([-0.5 0.5])
        ylim([-7 7])
        
        subplot(2,2,4)
        
        
        
        plot(currBinDataPts(:,3),mod(180/pi*atan2(currBinDataPts(:,2),currBinDataPts(:,1)),360),'.','Color',referenceColors(bi,:))
        %plot(distNorm,phasesDeg,'k.')
        hold on
        %dispVec=rotate2Dvector([currBinCenter currBinCenter],[pi/4]);
        
        %blueBar=plot([0 1],[0 360]-bi*distBinWidth,'b-','LineWidth',3)
        
        blueBar=plot([0 1],360-([angularBinCenters(bi) angularBinCenters(bi)]*180/pi),'b-','LineWidth',3);
       
        ylabel('Spike theta phase')
        xlabel('Distance in field')
        xlim([0 1])
        ylim([0 360])
        
        
        setFigFontTo(18)
        
        frameIdx=bi;
        
        try
            p1H.Visible='on';
            p2H.Visible='on';
            p3H.Visible='on';
            p4H.Visible='on';
            p5H.Visible='on';
        end
        
        try
            blueBar.Visible='on';
            
        end
        
         try
            pD.Visible='on';
        end
        
        F(frameIdx)=getframe(fCrossSections)
        
        try
            p1H.Visible='off';
            p2H.Visible='off';
            p3H.Visible='off';
            p4H.Visible='off';
            p5H.Visible='off';
        end
        
        try
            blueBar.Visible='off';
            
        end
        
         try
            pD.Visible='off';
        end
    end


    writerObj=VideoWriter(sprintf('%s.avi',titleStr));
    writerObj.FrameRate=1.5;
    
    open(writerObj);
    
    for i=1:length(F)
        try
        frame=F(i);
        
        writeVideo(writerObj,frame);
        end
    end
    close(writerObj)
    devFromLinearVsSpace=meanLinearDevPerBin;



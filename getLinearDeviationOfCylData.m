function [devFromLinearVsSpace] = getLinearDeviationOfCylData(phasesDegs,distNorm,distBinWidth)
%given a list of phases (0 to 360) and linear distances (between 0 and 1)
%returns the avg. deviation from linear relationship between phase and distance
%on 3D cylinder for each distance bin from 0 to 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %convert input into (x,y,z) triplets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phasesRadians=phasesDegs*pi/180;

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
    
    
    binCenters=(distBinWidth/2):distBinWidth:(1-distBinWidth/2);
    
    linearRefAngles=(1-binCenters)*2*pi; %linear phase precession points
    linearRefPts=[cos(linearRefAngles(:)),sin(linearRefAngles(:)),binCenters(:)];
    
    linearRefTangents=NaN(numBins,3);
    linearRefAxials=NaN(numBins,3);
    linearRefPerps=NaN(numBins,3);
    
    %figure
    for ri=1:(numBins)
        centerPt=[0 0 binCenters(ri)]
        if(ri+1>numBins)
            linearRefTangents(ri,:)=[0 0 0];
            linearRefAxials(ri,:)=[0 0 0];
        else
            linearRefTangents(ri,:)=linearRefPts(ri+1,:)-linearRefPts(ri,:);
            linearRefTangents(ri,:)=linearRefTangents(ri,:)/norm(linearRefTangents(ri,:));
            
             %towards center of circle is perpendicular, plus 1 radius higher 
            %gives perp vector in 45 degree angle plane for each z
            %linearRefAxials(ri,:)=(centerPt+[0 0 1])-linearRefPts(ri,:); 
            linearRefAxials(ri,:)=(centerPt)-linearRefPts(ri,:); 
            
            %linearRefAxials(ri,3)=0; %should be flat in z

            linearRefAxials(ri,:)=linearRefAxials(ri,:)/norm(linearRefAxials(ri,:));
        end
        %linearRefTangents(ri,:)=linearRefPts(ri+1,:);
        
        linearRefPerps(ri,:)=cross(linearRefTangents(ri,:),linearRefAxials(ri,:));
        
       
         %plot3Dvectors([linearRefTangents(ri,:);linearRefPerps(ri,:)]/10,linearRefPts(ri,:),'kr')
          %plot3DvectorFromTo(linearRefPts(ri,:),linearRefPts(ri+1,:),'k')
          %plot3DvectorFromTo(linearRefPts(ri,:),centerPt+[0 0 0.5],'b')
          
          %plot3DvectorFromTo(linearRefTangents(ri,:))
          %plot3DvectorFromTo(linearRefPerps(ri,:))
          
    end
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop through ellipitical bins
    %to store deviation from linear reference
    %for each distance bin
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    meanLinearDevPerBin=NaN(numBins,1);
    fCrossSections=figure;
    
    for i=1:numPts
        currPtPlot=pts(i,:);
        plot3(currPtPlot(1),currPtPlot(2),currPtPlot(3),'.')
        hold on
    end
    
    for bi=1:numBins
        currBinCenter=binCenters(bi); %same ellipse, only center shifting up
        currBinRefPt=linearRefPts(bi,:);
        
        %minorAxisUnitVector=[0 1 0];         %always parallel to base
        minorAxisUnitVector=[1 0 0];         %pooint to 0 degrees
         minorAxisUnitVector=currBinRefPt-[0 0 currBinCenter];         %pooint to 0 degrees
        %majorAxisUnitVector=[-1 0 1]/sqrt(2); %always 45degrees from base and 0 in minor direction
        %majorAxisUnitVector=[0 1 1]/sqrt(2); %always 45degrees from base and 0 in minor direction
        %majorAxisUnitVector=currBinRefPt-linearRefPerps(bi,:);
         majorAxisUnitVector=-linearRefPerps(bi,:);
        majorAxisUnitVector=majorAxisUnitVector/norm(majorAxisUnitVector);
        minorAxisUnitVector=minorAxisUnitVector/norm(minorAxisUnitVector);
        
        biFixed=1; %just use one major and minor axis, shifting only center vertically
        %biFixed=bi; 
        currBinUpperEllipseParams.centerPt=[0 0 currBinCenter+(distBinWidth/2)]; %same center as unit circle centered at z
        %currBinUpperEllipseParams.majorAxisUnitVector=linearRefPerps(biFixed,:);         %points perpendicular to ref. line
        currBinUpperEllipseParams.majorAxisUnitVector=majorAxisUnitVector;        %points perpendicular to ref. line
        %currBinUpperEllipseParams.minorAxisUnitVector=linearRefTangents(biFixed,:);      %points along to ref. line
        currBinUpperEllipseParams.minorAxisUnitVector=minorAxisUnitVector;      %points along to ref. line
        currBinUpperEllipseParams.majorAxisLength=sqrt(2);
        currBinUpperEllipseParams.minorAxisLength=1;
        currBinUpperEllipseParams.colorStr='k';
        
        currBinLowerEllipseParams.centerPt=[0 0 currBinCenter-(distBinWidth/2)]; %same center as unit circle centered at z
        %currBinLowerEllipseParams.majorAxisUnitVector=linearRefPerps(biFixed,:);         %points perpendicular to ref. line
        currBinLowerEllipseParams.majorAxisUnitVector=majorAxisUnitVector;       %points perpendicular to ref. line
        %currBinLowerEllipseParams.minorAxisUnitVector=linearRefTangents(biFixed,:);      %points along to ref. line
        currBinLowerEllipseParams.minorAxisUnitVector=minorAxisUnitVector;      %points along to ref. line
        currBinLowerEllipseParams.majorAxisLength=sqrt(2);
        currBinLowerEllipseParams.minorAxisLength=1;
        currBinLowerEllipseParams.colorStr='b';
        
        isAboveLowerEllipse=isAboveOrOnEllipse(pts,currBinLowerEllipseParams,fCrossSections);
        isAboveUpperEllipse=isAboveOrOnEllipse(pts,currBinUpperEllipseParams,fCrossSections);
        isBelowUpperEllipse=~isAboveUpperEllipse;
        
        isInCurrBin=isAboveLowerEllipse & isBelowUpperEllipse;
        
        
        
        currRefZ=currBinRefPt(3);
        
        currBinDataPts=pts(isInCurrBin,:);
        
        
        numPtsCurrBin=size(currBinDataPts,1);
        hold on
        scatter3(currBinDataPts(:,1),currBinDataPts(:,2),currBinDataPts(:,3),50*ones(numPtsCurrBin,1))
        scatter3(currBinRefPt(:,1),currBinRefPt(:,2),currBinRefPt(:,3),100)
        
        daspect([1 1 1])
        zlim([0 1])
        majorAxisVec=currBinLowerEllipseParams.majorAxisUnitVector*currBinLowerEllipseParams.majorAxisLength;
        minorAxisVec=currBinLowerEllipseParams.minorAxisUnitVector*currBinLowerEllipseParams.minorAxisLength;
        
        %plot3DvectorFromTo([0 0 0], linearRefPerps(biFixed,:),'b')
        %plot3DvectorFromTo([0 0 0], linearRefTangents(biFixed,:),'k')
        %plot3DvectorFromTo([0 0 0], majorAxisVec,'b')
        %plot3DvectorFromTo([0 0 0], minorAxisVec,'k')
        plot3DvectorFromTo(currBinRefPt, currBinRefPt-linearRefPerps(bi,:),'r')
        %plot3DvectorFromTo(currBinRefPt, currBinRefPt+linearRefAxials(bi,:),'r')
        
        
        
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
        
    end


    
    devFromLinearVsSpace=meanLinearDevPerBin;



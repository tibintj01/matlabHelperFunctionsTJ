function [] = plotPolarSpaceTimeField(spaceBins,timeBins,allFieldSpaceTimePhaseTriplets,fH)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    surfMode=0;
    polarMode=1;
    plotTimeAsColor=false;
    %plotTimeAsColor=true;
    %numCellsUsed=;
    minSampleCount=10;
    minSampleCount=1;
    
    
    
    minTime=2;
    maxTime=10;
      %maxTime=8;
      maxTime=15;
      
          numSpaceBins=50+1;
        numTimeBins=15+1;
        spaceBins=1:numSpaceBins;
        timeBins=1:numTimeBins;
        
    if(~exist('fH','var'))
        fH=figure;
    end
    
    if(surfMode)
        [X,Y]=meshgrid(spaceBins,timeBins);
        Z=avgPhasePerSpaceTimeBin';
        surf(X,Y,Z)
    elseif(polarMode==1)
        %roundPhases=round(avgPhasePerSpaceTimeBin);
        numPhaseBins=360;
        %numPhaseBins=120;
        numPhaseBins=20;
        numPhaseBins=24;
         %numPhaseBins=10;
        %numPhaseBins=30;
         numPhaseBins=40;
          %numPhaseBins=60;
        %numPhaseBins=15;
        %numPhaseBins=8;
        phaseEdges=linspace(0,360,numPhaseBins+1);
        phaseBinAssignments=discretize(allFieldSpaceTimePhaseTriplets(:,3),phaseEdges);
        
        %phaseAndSpaceToTime=NaN(length(spaceBins),numPhaseBins);
        timesAccumulationHappens=0;
        if(plotTimeAsColor)
            phaseAndSpaceToTime=NaN(length(spaceBins),numPhaseBins);
            phaseAndSpaceToTimeCount=zeros(length(spaceBins),numPhaseBins);
        else
            phaseAndSpaceToTime=NaN(length(timeBins),numPhaseBins);
            phaseAndSpaceToTimeCount=zeros(length(timeBins),numPhaseBins);
        end
        for phaseI=1:numPhaseBins
            spaceTimesPerPhase=find(phaseBinAssignments==phaseI);
            
            %associatedRows=NaN(length(spaceTimesPerPhase),1);
            %associatedCols=NaN(length(spaceTimesPerPhase),1);
            for i=1:length(spaceTimesPerPhase)
                %[row,col]=ind2sub(size(avgPhasePerSpaceTimeBin),spaceTimesPerPhase(i));
                spaceBin=allFieldSpaceTimePhaseTriplets(spaceTimesPerPhase(i),1)+1;
                timeBin=allFieldSpaceTimePhaseTriplets(spaceTimesPerPhase(i),2)+1;
                if(plotTimeAsColor)
                    [spaceBin timeBin]=swap(spaceBin,timeBin);
                else
                    if(timeBin<=minTime || timeBin>maxTime)
                        continue
                    end
                end
                %numCellsUsed=numCellsUsed+1;
                %if(isnan(phaseAndSpaceToTime(row,phaseI)))
                if(isnan(phaseAndSpaceToTime(timeBin,phaseI)))
                    %phaseAndSpaceToTime(row,phaseI)=col;
                    phaseAndSpaceToTime(timeBin,phaseI)=spaceBin;
                else
                    %phaseAndSpaceToTime(row,phaseI)=nanmean([phaseAndSpaceToTime(row,phaseI); col]);
                    phaseAndSpaceToTime(timeBin,phaseI)=phaseAndSpaceToTime(timeBin,phaseI)+spaceBin;
                    %if(timeBin==14)
                     %   disp('here')  
                    %end
                    %timesAccumulationHappens=timesAccumulationHappens+1;
                end
                phaseAndSpaceToTimeCount(timeBin,phaseI)=phaseAndSpaceToTimeCount(timeBin,phaseI)+1;
                %associatedRows(i)=row;
                %associatedCols(i)=col;
            end
            %meanAssociatedRow=nanmean(associatedRows);
            %meanAssociatedCol=nanmean(associatedCols);
            %phaseAndSpaceToTime(round(meanAssociatedCol),phaseI)=meanAssociatedRow;
            
        end
       
        phaseAndSpaceToTimeCount(phaseAndSpaceToTimeCount==0)=NaN;
         lowCountSpaceTimePts=phaseAndSpaceToTimeCount<minSampleCount;
         
        phaseAndSpaceToTime=phaseAndSpaceToTime./phaseAndSpaceToTimeCount;
        phaseAndSpaceToTime(lowCountSpaceTimePts)=NaN;
        scaledPosition=scaledata(phaseAndSpaceToTime,0,1);
        %polarplot3d(scaledPosition,'PlotType','surf')
        polarplot3d(scaledPosition)
        %polarplot3d(scaledPosition,'PlotType','contour','ContourLines',500)
        colormap(jet)
        %colormap(parula)
        timesAccumulationHappens
        minDisp=prctile(scaledPosition(:),5);
        maxDisp=prctile(scaledPosition(:),99);
        %caxis([minDisp maxDisp])
        caxis([0 1])
        %title(sprintf('N=%d cells',numCellsUsed))
        
    end
    
   
    
        
    %cmp=hsv
    %hsvInvert=[cmp(ceil(end/2):end,:); cmp(1:floor(end/2),:)];
    %colormap(hsv)
    %colormap(hsvInvert)
    %cb=colorbar('north')
    %cb.Position=[0.05 0.25 0.4 0.05 ];
    cb=colorbar
    %xlabel('Pos (frac of track length)')
    %ylabel('Time since leaving wall (sec)')
    if(polarMode)
        if(plotTimeAsColor)
          ylabel(cb,'Time elapsed in field')
        else
          ylabel(cb,'Position (frac. of field)')
        end
    else
        ylabel(cb,'First spike theta phase (deg)')
    end
    %xlim([0 1])
    %ylim([-0.5 max(timeBins)+0.5])
    daspect([1 1 1])
    view([90 90])
    %{
    figure
    phaseBins=edgesToBins(phaseEdges);
      %   omarPcolor(phaseBins,timeBins,scaledPosition,fH)
       % colormap(jet)
        %figure
        repeatedPhaseBins=[phaseBins phaseBins+360];
        [X,Y]=meshgrid(repeatedPhaseBins,timeBins);
        %Z=avgPhasePerSpaceTimeBin';
        repeatedPosition=[scaledPosition,scaledPosition];
        %surf(X,Y,scaledPosition)
        K=1/9*(ones(3));
        
        smoothedRepeatedPos=conv2(repeatedPosition,K,'same');
        smoothedRepeatedPos=conv2(smoothedRepeatedPos,K,'same');
        smoothedRepeatedPos=conv2(smoothedRepeatedPos,K,'same');
        %surf(X,Y,repeatedPosition)
        subplot(2,1,1)
        surf(X,Y,smoothedRepeatedPos)
        xlim([160 160+310])
        subplot(2,1,2)
        contour(X,Y,smoothedRepeatedPos,25)
        %contour(X,Y,smoothedRepeatedPos,15)

        %xlim([150 150+360])
        xlim([150 150+310])
                disp('here')
    %}
                
end


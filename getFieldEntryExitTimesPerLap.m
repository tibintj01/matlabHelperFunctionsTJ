function [numLaps,numFields,fieldEntryTimes,fieldExitTimes] = getFieldEntryExitTimesPerLap(unitDataFilePath)
	
    unitData=load(unitDataFilePath);

    %showPlots=1;
        showPlots=1;
        showLapPlots=1;
        
        showPlots=0;
        showLapPlots=0;
        
        entryTimeIdxesPerDirPerFieldPerLap=[];
        exitTimeIdxesPerDirPerFieldPerLap=[];

        entryTimePerDirPerFieldPerLap=[];
        exitTimePerDirPerFieldPerLap=[];

        entryPosPerDirPerFieldPerLap=[];
        exitPosPerDirPerFieldPerLap=[];

        hasLapFieldInterference=[];
        
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%load relevant variables for computing field exit and entry times
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	positionTimeAxis=unitData.positionTimeAxis;
    approxPosTimeStep=median(diff(positionTimeAxis));
	positionPerTime=unitData.positionPerTimeStep;
    
    approxTimeStep=median(diff(positionTimeAxis));
    
    nonNanPosIdxes=find(~isnan(positionPerTime));
    firstNonNanPosIdx=nonNanPosIdxes(1);
   	lastNonNanPosIdx=nonNanPosIdxes(end);
    

	signedSpeedPerTime=unitData.speedPerTimeStep;

	hasLeftwardFields=0;
	hasRightwardFields=0;

	if(isfield(unitData,'manualFieldEndsMleftward'))
		manualFieldEndsMleftward=unitData.manualFieldEndsMleftward;
		hasLeftwardFields=1;
	end
	if(isfield(unitData,'manualFieldStartsMleftward'))
		manualFieldStartsMleftward=unitData.manualFieldStartsMleftward;
		hasLeftwardFields=1;
	end
	if(isfield(unitData,'manualFieldEndsMrightward'))
		manualFieldEndsMrightward=unitData.manualFieldEndsMrightward;
		hasRightwardFields=1;
	end
	if(isfield(unitData,'manualFieldStartsMrightward'))
		manualFieldStartsMrightward=unitData.manualFieldStartsMrightward;
		hasRightwardFields=1;
    end


    numLeftFields=0;
    numRightFields=0;
    if(hasLeftwardFields)
        numLeftFields=length(manualFieldStartsMleftward);
    end
    if(hasRightwardFields)
        numRightFields=length(manualFieldStartsMrightward);
    end
	
	%posClosenessThresh=0.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%for circular maze, circ shift position so that lowest speed position is 0
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isCircMaze=0;
	if(strcmp('Circular Maze',unitData.unitInfo.positionPerTime.MazeType))
		medianPos=nanmedian(unitData.positionPerTimeStep);
        %circumference=max(unitData.positionPerTimeStep);
        circumference=max(unitData.positionPerTimeStep((end/4):(3*end/4))); %ends are weird after filling
        isCircMaze=1;
	%circShiftPos=mod(unitData.positionPerTimeStep-medianPos,max(unitData.positionPerTimeStep));
		%circShiftPos=mod(unitData.positionPerTimeStep+medianPos,max(unitData.positionPerTimeStep));
    circShiftPos=unitData.positionPerTimeStep;
        
        positionPerTime=circShiftPos;	
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %recover lost position phase information by shifting until maximal
        %phase correlation in first field bounds
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %{
        original1DposPerTime=unitData.unitInfo.positionPerTime.OneDLocation;
        bestCircShift=0;
        badIdxes=isnan(original1DposPerTime);
        zeroedOriginal=zscoreLFP(original1DposPerTime);
        %}
        
        savedPosPerTimeStep=unitData.positionPerTimeStep;
        %circ always leftward!
        spikePhasePerTimeStep=unitData.directionSpecificStats.leftward.localSpikePhasePerTime;
        
         if(~contains(unitDataFilePath,'Gatsby'))
                spikePhasePerTimeStep=mod(spikePhasePerTimeStep+180,360);
          end
        
        firstFieldStartSavedPos=unitData.manualFieldStartsMleftward(1);
        firstFieldEndSavedPos=unitData.manualFieldEndsMleftward(1);
        
        leftSpikePositions=unitData.leftSpikePositions;
        
        numShifts=500;
        thisShiftCorr=NaN(numShifts,1);
        for csi=1:numShifts
            currShift=(csi/numShifts)*circumference;
            currShiftedPos=mod(savedPosPerTimeStep+currShift,circumference);
            invalidIdxes=currShiftedPos>firstFieldStartSavedPos | currShiftedPos<firstFieldEndSavedPos;
            currShiftedPos(invalidIdxes)=NaN;
            currPhases=spikePhasePerTimeStep;
            currPhases(invalidIdxes)=NaN;
            
            inFieldSpikeCount=sum(~isnan(currPhases));
            
            %if(sum(~isnan(currPhases))<3000)
            %    continue
            %end
            
            thisShiftCorr(csi)=nancorr(currShiftedPos,currPhases)*inFieldSpikeCount;
            %thisShiftCorr(csi)=inFieldSpikeCount;
         
           
            %if(csi==100)
            %   figure;
            %    plot(currShiftedPos,currPhases,'k.') 
            %end
            %matchForShift(csi)=nansum(zeroedOriginal.*zscoreLFP(savedVersion));
            %matchForShift(csi)=nancorr(zeroedOriginal,zscoreLFP(shiftedPos));
             %matchForShift(csi)=nansum(abs(original1DposPerTime-shiftedPos));
        end
          [~,bestCircShiftId]=max((thisShiftCorr));
          bestCircShift=(bestCircShiftId/numShifts)*circumference;
        
        %savedVersion(badIdxes)=NaN;
        %savedVersion=zscoreLFP(savedVersion);
        
        %originalCircum=max(original1DposPerTime);
        %savedVersion=savedVersion/max(savedVersion)*max(original1DposPerTime);
        %{
        numShifts=1000;
        for csi=1:numShifts
            currShift=(csi/numShifts)*circumference;
            shiftedPos=mod(savedVersion+currShift,circumference);
            %matchForShift(csi)=nansum(zeroedOriginal.*zscoreLFP(savedVersion));
            %matchForShift(csi)=nancorr(zeroedOriginal,zscoreLFP(shiftedPos));
             matchForShift(csi)=nansum(abs(original1DposPerTime-shiftedPos));
        end
          [~,bestCircShiftId]=min(matchForShift);
          bestCircShift=(bestCircShiftId/numShifts)*3;
          
          %}
          
          circShiftPos=mod(unitData.positionPerTimeStep+bestCircShift,circumference);
            positionPerTime=circShiftPos;
            
            circShiftPos=fillmissing(circShiftPos,'linear');
            positionPerTime=fillmissing(positionPerTime,'linear');
         
         %
         %{
    figure; plot(positionPerTime,spikePhasePerTimeStep,'k.')
        hold on
        plot([firstFieldStartSavedPos firstFieldStartSavedPos],ylim,'g--')
        plot([firstFieldEndSavedPos firstFieldEndSavedPos],ylim,'r--')
        %plot(positionTimeAxis,positionPerTime,'r-','LineWidth',2)
        figure; plot(thisShiftCorr)
        disp('')
         %}
        
            
        
        
            
        
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %shift entered field boundaries accordingly
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(hasRightwardFields)
            manualFieldStartsMrightward=mod(manualFieldStartsMrightward-medianPos,max(unitData.positionPerTimeStep));
            manualFieldEndsMrightward=mod(manualFieldEndsMrightward-medianPos,max(unitData.positionPerTimeStep));
        elseif(hasLeftwardFields)
            manualFieldStartsMleftward=mod(manualFieldStartsMleftward-medianPos,max(unitData.positionPerTimeStep));
            manualFieldEndsMleftward=mod(manualFieldEndsMleftward-medianPos,max(unitData.positionPerTimeStep));
        end
        %}
        
		%positionPerTime=scaledata(positionPerTime,0,pi); %2.9 was because of cutting of ends in provided 1d data

		approxPosTimeStep=median(diff(positionTimeAxis));

			polynomialFitOrder=1; %fit lines over acc time width
		    speedEstTimeWidth=0.125; %sec %~5 point linear fit estimation of speed
		      %speedEstTimeWidth=0.075; %sec %~5 point linear fit estimation of speed
		    supportWindLength=ceil(speedEstTimeWidth/approxPosTimeStep);
            %unwrap considering position as radians (needs to be times 2, divide out later)
            unwrappedSpeedPerTimeStep=movingslope(unwrap(circShiftPos*2),supportWindLength,polynomialFitOrder,approxPosTimeStep);
            unwrappedSpeedPerTimeStep=unwrappedSpeedPerTimeStep/2;
            
		    recalcAbsSpeedPerTimeStep=abs(unwrappedSpeedPerTimeStep); %meters/sec

        %{
		[~,outlierIdxes,~] = deleteoutliers(recalcSpeedPerTimeStep);
		speedWithoutCircTeleports=recalcSpeedPerTimeStep;
		speedWithoutCircTeleports(outlierIdxes)=NaN;
        signedSpeedPerTime=speedWithoutCircTeleports;
        %}
       signedSpeedPerTime=unwrappedSpeedPerTimeStep;
        %filledPositionPerTime=positionPerTime;

    else
        %filledPositionPerTime=fillmissing(positionPerTime,'linear');
        %filledPositionPerTime=positionPerTime;
        unwrappedSpeedPerTimeStep=signedSpeedPerTime;
        medianPos=NaN;
    end
	
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find large regions of consecutive NaN and label so that they aren't
    %filled
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nanLocations = isnan(positionPerTime) % Get logical array of whether element is NaN or not.
    props = regionprops(nanLocations, 'Area', 'PixelIdxList'); % Find all the regions.
    
    timeThreshForFilling=6; %sec
     timeThreshForFilling=10; %sec
    doNotFillSentinelVal=-100;
    
    for k = 1 : length(props)
        currentBlockDuration=props(k).Area*approxTimeStep;
        currentBlockStartIdx=props(k).PixelIdxList(1);
        currentBlockEndIdx=props(k).PixelIdxList(end);
  fprintf('Region #%d has length %d sec and starts at element %d and ends at element %d\n',...
    k, currentBlockDuration, currentBlockStartIdx, currentBlockEndIdx);

        if(currentBlockDuration>timeThreshForFilling)
            positionPerTime(currentBlockStartIdx:currentBlockEndIdx)=doNotFillSentinelVal;
        end
    end
    

    filledPositionPerTime=fillmissing(positionPerTime,'linear');
    filledPositionPerTime(filledPositionPerTime==doNotFillSentinelVal)=NaN;
    positionPerTime(positionPerTime==doNotFillSentinelVal)=NaN;
    

			polynomialFitOrder=1; %fit lines over acc time width
		    speedEstTimeWidth=0.125; %sec %~5 point linear fit estimation of speed
		      %speedEstTimeWidth=0.075; %sec %~5 point linear fit estimation of speed
		    supportWindLength=ceil(speedEstTimeWidth/approxPosTimeStep);
            %unwrap considering position as radians (needs to be times 2, divide out later)
            signedSpeedPerTime=movingslope(filledPositionPerTime,supportWindLength,polynomialFitOrder,approxPosTimeStep);
    
    %figure; plot(positionTimeAxis,filledPositionPerTime)
        
    leftwardPositionPerTime=positionPerTime; leftwardPositionPerTime(signedSpeedPerTime>0)=NaN;
	rightwardPositionPerTime=positionPerTime; rightwardPositionPerTime(signedSpeedPerTime<0)=NaN;
    
    leftwardSpeedIdxes=signedSpeedPerTime<0;
    rightwardSpeedIdxes=signedSpeedPerTime>0;
	
    filledPositionPerTime(1:(firstNonNanPosIdx-1))=NaN;
    filledPositionPerTime((lastNonNanPosIdx+1):end)=NaN;

	%filledSignedSpeedPerTime=fillmissing(signedSpeedPerTime,'linear');
    filledSignedSpeedPerTime=signedSpeedPerTime;
    
    filledSignedSpeedPerTime(1:(firstNonNanPosIdx-1))=NaN;
    filledSignedSpeedPerTime((lastNonNanPosIdx+1):end)=NaN;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%segment position data (in current direction) into laps by starting in middle and searching outwards
	%for linear maze, lap segmented by first time rat turns around past 95% of track length
	%for circular maze, laps segmented by first time the position jumps a full track within a time step
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    closenessThresh=0.02; %2 cm
    
    
    middlePosition=min(positionPerTime)+range(filledPositionPerTime)/2;
    lowerPositions=positionPerTime;
    lowerPositions(positionPerTime>middlePosition)=NaN;
    
    upperPositions=positionPerTime;
    upperPositions(positionPerTime<middlePosition)=NaN;
    
    [lowerPosCounts,~,discLowerPos]=histcounts(lowerPositions,100);
    [upperPosCounts,~,discUpperPos]=histcounts(upperPositions,100);
    
    [~,statisticalPosMinBin]=max(lowerPosCounts);
    [~,statisticalPosMaxBin]=max(upperPosCounts);
    
    if(isCircMaze)
        statisticalPosMax=circumference*0.95;
        statisticalPosMin=circumference*0.05;
    else
        statisticalPosMin=nanmean(positionPerTime(discLowerPos==statisticalPosMinBin));
        statisticalPosMax=nanmean(positionPerTime(discUpperPos==statisticalPosMaxBin));
    end
    
    statisticalRange=statisticalPosMax-statisticalPosMin;
    statisticalMiddlePos=statisticalPosMin+statisticalRange/2;
    
    lapSeedPosIdxes=find(abs(filledPositionPerTime-statisticalMiddlePos)<closenessThresh);
       %maxPositionHigh=statisticalPosMin+0.95*statisticalRange;
        %maxPositionLow=statisticalPosMin+0.05*statisticalRange;
        
        maxPositionHigh=statisticalPosMin+0.975*statisticalRange;
        maxPositionLow=statisticalPosMin+0.025*statisticalRange;
        
        %maxPositionHigh=statisticalPosMin+0.925*statisticalRange;
        %maxPositionLow=statisticalPosMin+0.075*statisticalRange;
        
    if(showPlots)
        close all
        figure; plot(positionTimeAxis,filledPositionPerTime); hold on
        plot(positionTimeAxis,filledSignedSpeedPerTime)
          %figure; plot(positionTimeAxis,positionPerTime); hold on
        %plot(positionTimeAxis,signedSpeedPerTime)
        plot(positionTimeAxis(lapSeedPosIdxes),middlePosition,'ko','MarkerSize',5)

        plot(xlim,[maxPositionHigh maxPositionHigh],'k--')
        plot(xlim,[maxPositionLow maxPositionLow],'k--')
        

            plot(xlim,[0 0],'r--')
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %leftward movement case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    leftwardLapStartIdxesNonUnique=[];
    leftwardLapEndIdxesNonUnique=[];
    rightwardLapStartIdxesNonUnique=[];
    rightwardLapEndIdxesNonUnique=[];
    
    for si=1:length(lapSeedPosIdxes)
        currLapSeedPosIdx=lapSeedPosIdxes(si);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find lap start relative to this seed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        currPos=filledPositionPerTime(currLapSeedPosIdx);
        currSpeed=filledSignedSpeedPerTime(currLapSeedPosIdx);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %determine if this is rightward or leftward lap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        findLapDirectionOfCurrSeed

        if(currSeedInLeftwardLap)
            setCurrLapBoundsLeftwardLap
            leftwardLapStartIdxesNonUnique=[leftwardLapStartIdxesNonUnique; firstIdxOfThisLap];
            leftwardLapEndIdxesNonUnique=[leftwardLapEndIdxesNonUnique; lastIdxOfThisLap];
        else
            setCurrLapBoundsRightwardLap
            rightwardLapStartIdxesNonUnique=[rightwardLapStartIdxesNonUnique; firstIdxOfThisLap];
            rightwardLapEndIdxesNonUnique=[rightwardLapEndIdxesNonUnique; lastIdxOfThisLap];
        end
        
        
    end
    
   
    leftwardLapStartStopIdxes=unique([leftwardLapStartIdxesNonUnique leftwardLapEndIdxesNonUnique],'rows');

    rightwardLapStartStopIdxes=unique([rightwardLapStartIdxesNonUnique rightwardLapEndIdxesNonUnique],'rows');
    
    removeRows=[];
    for li=2:size(leftwardLapStartStopIdxes,1)
        if(size(leftwardLapStartStopIdxes,1) >1 &&  leftwardLapStartStopIdxes(li,1)==leftwardLapStartStopIdxes(li-1,1) || leftwardLapStartStopIdxes(li,2)==leftwardLapStartStopIdxes(li-1,2))
            removeRows=[removeRows; li];
        end
    end
    leftwardLapStartStopIdxes(removeRows,:)=[];
    
    removeRows=[];
    for li=2:size(rightwardLapStartStopIdxes,1)
        if(size(rightwardLapStartStopIdxes,1) >1 && rightwardLapStartStopIdxes(li,1)==rightwardLapStartStopIdxes(li-1,1) || rightwardLapStartStopIdxes(li,2)==rightwardLapStartStopIdxes(li-1,2))
            removeRows=[removeRows; li];
        end
    end
    rightwardLapStartStopIdxes(removeRows,:)=[];

    if(showPlots)
       for li=1:length(leftwardLapStartStopIdxes)

           currPosTime=positionTimeAxis(leftwardLapStartStopIdxes(li,1));
           plot([currPosTime currPosTime],ylim,'b-','LineWidth',4)


            currPosTime=positionTimeAxis(leftwardLapStartStopIdxes(li,2));
           plot([currPosTime currPosTime],ylim,'k-','LineWidth',4)

       end
    end
	
    
    if(~isempty(leftwardLapStartStopIdxes))
        leftwardLapDurations=(leftwardLapStartStopIdxes(:,2)-leftwardLapStartStopIdxes(:,1))*approxPosTimeStep;
        fakeLapIDs=leftwardLapDurations<1; %1 second is unreasonably fast for a whole lap
        leftwardLapStartStopIdxes(fakeLapIDs,:)=[];
        leftwardLapDurations(fakeLapIDs)=[];
         if(~isempty(leftwardLapDurations) && ~isempty(deleteoutliers(leftwardLapDurations)))
            [~,outlierDurIdxLeft,~]=deleteoutliers(leftwardLapDurations);
              leftwardLapStartStopIdxes(outlierDurIdxLeft,:)=[];
                  leftwardLapDurations(outlierDurIdxLeft)=[];
         end
                numLapsLeftward=length(leftwardLapDurations);
    else
         numLapsLeftward=0;
     end
        
     if(~isempty(rightwardLapStartStopIdxes))
        rightwardLapDurations=(rightwardLapStartStopIdxes(:,2)-rightwardLapStartStopIdxes(:,1))*approxPosTimeStep;
        fakeLapIDs=rightwardLapDurations<1; %1 second is unreasonably fast for a whole lap
        rightwardLapStartStopIdxes(fakeLapIDs,:)=[];
         rightwardLapDurations(fakeLapIDs)=[];

         if(~isempty(rightwardLapDurations) && ~isempty(deleteoutliers(rightwardLapDurations)))
       [~,outlierDurIdxRight,~]=deleteoutliers(rightwardLapDurations);
        rightwardLapStartStopIdxes(outlierDurIdxRight,:)=[];
        rightwardLapDurations(outlierDurIdxRight)=[];
         end



        numLapsRightwards=length(rightwardLapDurations);
      else
         numLapsRightwards=0;
     end
    
  disp('')
    
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %loop through both directions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        maxNumFields=5;
        maxNumLaps=1000;
        %entryTimeIdxesPerDirPerFieldPerLap=NaN(2,maxNumFields,maxNumLaps);
        %exitTimeIdxesPerDirPerFieldPerLap=NaN(2,maxNumFields,maxNumLaps);

        %entryTimePerDirPerFieldPerLap=NaN(2,maxNumFields,maxNumLaps);
       
        %fieldExitTimesPerLap=NaN(2,maxNumFields,maxNumLaps);
        numFieldsPerDir=NaN(2,1);
        numLapsPerDir=NaN(2,1);
        
        
        if(showLapPlots)
            figure; 
        end
         %hasLapFieldInterference=zeros(2,numFields,currNumLaps);
         
        for di=1:2
            if(di==1)
                if(~hasRightwardFields)
                    continue
                end
                numFields=numRightFields;
                currFieldStarts=manualFieldStartsMrightward;
                currFieldEnds=manualFieldEndsMrightward;
                currPosPerTime=rightwardPositionPerTime;
                currNumLaps=numLapsRightwards;
                currDirLapStartStopIdxes=rightwardLapStartStopIdxes;
                currDirLapStartTimes=positionTimeAxis(rightwardLapStartStopIdxes(:,1));
                currDirLapStopTimes=positionTimeAxis(rightwardLapStartStopIdxes(:,2));
                dirStr='rightward';
            elseif(di==2)
                if(~hasLeftwardFields)
                    continue
                end
                numFields=numLeftFields;
                currFieldStarts=manualFieldStartsMleftward;
                currFieldEnds=manualFieldEndsMleftward;
                currPosPerTime=leftwardPositionPerTime;
                currNumLaps=numLapsLeftward;
                 currDirLapStartStopIdxes=leftwardLapStartStopIdxes;
                  currDirLapStartTimes=positionTimeAxis(leftwardLapStartStopIdxes(:,1));
                currDirLapStopTimes=positionTimeAxis(leftwardLapStartStopIdxes(:,2));
                 dirStr='leftward';
            end
            
            numFieldsPerDir(di)=numFields;
            numLapsPerDir(di)=currNumLaps;
            lapStartTimesPerDir.(dirStr)=currDirLapStartTimes;
            lapStopTimesPerDir.(dirStr)=currDirLapStopTimes;            


	     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	     %adjust each lap start end time  to include when field lies beyond lap boundary
	     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %figure
            for fi=1:numFields
                currFieldStartPos=currFieldStarts(fi);
                currFieldEndPos=currFieldEnds(fi);
            	for li=1:currNumLaps
		    currLapStartTimeIdx=currDirLapStartStopIdxes(li,1);
                    currLapEndTimeIdx=currDirLapStartStopIdxes(li,2);		
	    	    thisLapPositions=currPosPerTime(currLapStartTimeIdx:currLapEndTimeIdx);
                
                thisLapPositionsNoNan=thisLapPositions(~isnan(thisLapPositions));
                    if(isempty(thisLapPositionsNoNan))
                        continue
                    end
                    
                    currLapStartPos=thisLapPositionsNoNan(1);
                    currLapEndPos=thisLapPositionsNoNan(end);
                    
	
		     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		     %if field bound is 'outside' lap bound, extend lap bound to field bound
		     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        lapStartMinusFieldStart=currLapStartPos-currFieldStartPos;
                        lapEndMinusFieldEnd=currLapEndPos-currFieldEndPos;
		    fieldStartsBeforeLap=false;
		    fieldEndsAfterLap=false;
		    extendLapBackwardTimeIdx=currLapStartTimeIdx;
		    extendLapFwdTimeIdx=currLapEndTimeIdx;
            fieldTolerance=0.01;
		    if(di==1) %rightward, increasing position
			if(lapStartMinusFieldStart>0)
				fieldStartsBeforeLap=true;
				extendLapBackwardTimeIdx=currLapStartTimeIdx;
				while(currPosPerTime(extendLapBackwardTimeIdx)- currFieldStartPos >fieldTolerance)
					extendLapBackwardTimeIdx=extendLapBackwardTimeIdx-1;
					if(extendLapBackwardTimeIdx==1)
						break;
					end
				end
			end
			if(lapEndMinusFieldEnd<0)
                                fieldEndsAfterLap=true;
				extendLapFwdTimeIdx=currLapEndTimeIdx;
                                while(currPosPerTime(extendLapFwdTimeIdx)- currFieldEndPos <-fieldTolerance)
                                        extendLapFwdTimeIdx=extendLapFwdTimeIdx+1;
                                        if(extendLapFwdTimeIdx==length(currPosPerTime))
                                                break;
                                        end
                                end
                        end
           elseif(di==2)%leftward, decreasing position
			if(lapStartMinusFieldStart<0)
				fieldStartsBeforeLap=true;
                                extendLapBackwardTimeIdx=currLapStartTimeIdx;
                                while(currPosPerTime(extendLapBackwardTimeIdx)- currFieldStartPos <-fieldTolerance)
                                        extendLapBackwardTimeIdx=extendLapBackwardTimeIdx-1;
                                        if(extendLapBackwardTimeIdx==1)
                                                break;
                                        end
                                end
			end
			if(lapEndMinusFieldEnd>0)
                                fieldEndsAfterLap=true;
                                extendLapFwdTimeIdx=currLapEndTimeIdx;
                                while( currPosPerTime(extendLapFwdTimeIdx)- currFieldEndPos >fieldTolerance)
                                        extendLapFwdTimeIdx=extendLapFwdTimeIdx+1;
                                        if(extendLapFwdTimeIdx==length(currPosPerTime))
                                                break;
                                        end
                                end
                        end
                    end

		    
                %currDirLapStartStopIdxes(li,1)=extendLapBackwardTimeIdx;
		    %currDirLapStartStopIdxes(li,2)=extendLapFwdTimeIdx;
            thisLapPositions=currPosPerTime( currDirLapStartStopIdxes(li,1):currDirLapStartStopIdxes(li,2));
            %plot(thisLapPositions,'ko')
            %hold on
		end	
	    end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %loop through each available field boundary set
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            for fi=1:numFields
                
                if(isCircMaze)
                    %currStartPos=mod(currFieldStarts(fi)-medianPos,circumference);
                    %currEndPos=mod(currFieldEnds(fi)-medianPos,circumference);
                    currStartPos=currFieldStarts(fi);
                    currEndPos=currFieldEnds(fi);
                else
                    currStartPos=currFieldStarts(fi);
                    currEndPos=currFieldEnds(fi);
                end
  
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %loop through each lap
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for li=1:currNumLaps
                    currLapStartTimeIdx=currDirLapStartStopIdxes(li,1);
                    currLapEndTimeIdx=currDirLapStartStopIdxes(li,2);
                    
                    thisLapPositions=currPosPerTime(currLapStartTimeIdx:currLapEndTimeIdx);

                    thisLapPositionsNoNan=thisLapPositions(~isnan(thisLapPositions));
                    if(isempty(thisLapPositionsNoNan))
                        continue
                    end
                    
                    currLapStartPos=thisLapPositionsNoNan(1);
                    currLapEndPos=thisLapPositionsNoNan(end);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %lap end takes precedence over indicated field end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if(di==1) %rightward, increasing position
                        effectiveFieldStart=max(currLapStartPos,currStartPos);
                        effectiveFieldEnd=min(currLapEndPos,currEndPos);
                    elseif(di==2)%leftward, decreasing position
                        effectiveFieldStart=min(currLapStartPos,currStartPos);
                        effectiveFieldEnd=max(currLapEndPos,currEndPos);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %already adjusted laps
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    effectiveFieldStart=currStartPos;
                    effectiveFieldEnd=currEndPos;
                    
                    if(isCircMaze)
                        thisLapPositionsUnwrapped=unwrap(thisLapPositions*2)/2;
                    else
                        thisLapPositionsUnwrapped=thisLapPositions;
                    end
                    
                    %okay to fill within lap?
                    %thisLapPositionsFilledAndUnwrapped=fillmissing(thisLapPositionsUnwrapped,'linear');
                    thisLapPositionsFilledAndUnwrapped=thisLapPositionsUnwrapped;
                    thisLapPositionsFilledAndUnwrapped=scaledata(thisLapPositionsFilledAndUnwrapped,min(thisLapPositions),max(thisLapPositions));
                    if(showLapPlots)
                        hold on
                        plot(positionTimeAxis(currLapStartTimeIdx:currLapEndTimeIdx),thisLapPositionsFilledAndUnwrapped,'ko')
                        title(unitData.unitInfo.sessionName)
                        %shownTitles={shownTitles, unitData.unitInfo.sessionName};
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %find initial times corresponding to start and end
                    %crossings CIRC DIFF FOR CIRC MAZE
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %{
                        if(unitData.isCircMaze)
                        effectiveFieldStartArray=repmat(effectiveFieldStart,size(thisLapPositionsFilledAndUnwrapped));
                        effectiveFieldEndArray=repmat(effectiveFieldEnd,size(thisLapPositionsFilledAndUnwrapped));
                        isCircCloseStart=abs(circDiffGivenCircumference(thisLapPositionsFilledAndUnwrapped,effectiveFieldStartArray,max(unitData.positionPerTimeStep)))<=closenessThresh;
                        isCircCloseEnd=abs(circDiffGivenCircumference(thisLapPositionsFilledAndUnwrapped,effectiveFieldEndArray,max(unitData.positionPerTimeStep)))<=closenessThresh;
                        closeEntryIdxes=find(isCircCloseStart);
                        closeExitIdxes=find(isCircCloseEnd);
                        else
                     %}
                        closeEntryIdxes=find(abs(thisLapPositionsFilledAndUnwrapped-effectiveFieldStart)<=closenessThresh);
                        closeExitIdxes=find(abs(thisLapPositionsFilledAndUnwrapped-effectiveFieldEnd)<=closenessThresh);
                    %end
                    
                    if(isempty(closeEntryIdxes))
                        hasLapFieldInterference(di,fi,li)=1;
                        %fieldStartIdxTemp=1;
                        continue
                    else
                        
                        fieldStartIdxTemp=min(closeEntryIdxes);
                    end
                     if(isempty(closeExitIdxes))
                       %fieldEndIdxTemp=currLapEndTimeIdx;
                       hasLapFieldInterference(di,fi,li)=1;
                       continue
                     else
                         fieldEndIdxTemp=max(closeExitIdxes);
                     end
                    hasLapFieldInterference(di,fi,li)=0;
                  
                    %currFieldFirstEntryTimeIdx=min(closeEntryIdxes)+currLapStartTimeIdx-1;
                    %currFieldLastExitTimeIdx=max(closeExitIdxes)+currLapStartTimeIdx-1;
                      currFieldFirstEntryTimeIdx=fieldStartIdxTemp+currLapStartTimeIdx-1;
                    currFieldLastExitTimeIdx=fieldEndIdxTemp+currLapStartTimeIdx-1;
                    
                    if(currFieldFirstEntryTimeIdx<1)
                        currFieldFirstEntryTimeIdx=1;
                    end
                    
                    if(currFieldLastExitTimeIdx>length(positionTimeAxis))
                        currFieldLastExitTimeIdx=length(positionTimeAxis);
                    end
                    
                    entryTimeIdxesPerDirPerFieldPerLap(di,fi,li)=currFieldFirstEntryTimeIdx;
                    exitTimeIdxesPerDirPerFieldPerLap(di,fi,li)=currFieldLastExitTimeIdx;
                    entryTimePerDirPerFieldPerLap(di,fi,li)=positionTimeAxis(currFieldFirstEntryTimeIdx);
                    exitTimePerDirPerFieldPerLap(di,fi,li)=positionTimeAxis(currFieldLastExitTimeIdx);
                    
                    entryPosPerDirPerFieldPerLap(di,fi,li)=positionPerTime(currFieldFirstEntryTimeIdx);
                    exitPosPerDirPerFieldPerLap(di,fi,li)=positionPerTime(currFieldLastExitTimeIdx);
                    
                    if(showLapPlots)
                        currEntryPos=positionPerTime(currFieldFirstEntryTimeIdx);
                        currExitPos=positionPerTime(currFieldLastExitTimeIdx);
                        currEntryTime=positionTimeAxis(currFieldFirstEntryTimeIdx);
                        currExitTime=positionTimeAxis(currFieldLastExitTimeIdx);
                        hold on
                        
                        plot(positionTimeAxis(currLapStartTimeIdx:currLapEndTimeIdx), positionPerTime(currLapStartTimeIdx:currLapEndTimeIdx))
                        
                        fieldWidths(fi,li)=abs(currExitPos-currEntryPos);
                        
                        %if(abs(currExitPos-currEntryPos)<0.2)
                        %    disp('')   
                        %end
                        
                        plot(xlim,[currEntryPos currEntryPos],'g--','LineWidth',3)
                        plot(xlim,[currExitPos currExitPos],'r--','LineWidth',3)
                        
                         plot([currEntryTime currEntryTime],ylim,'g--','LineWidth',3)
                        plot([currExitTime currExitTime],ylim,'r--','LineWidth',3)
                        
                        
                        %plot(positionTimeAxis,positionPerTime,'mo')
                        
                        %plotOriginalSpikeTimes(unitData)
                                                
                                                    
                     
                    end
                    
                end %lap loop
                %{
                   for s=1:length(unitData.leftSpikeTimes(:))
                             plot([unitData.leftSpikeTimes(s) unitData.leftSpikeTimes(s)],ylim,'k-','LineWidth',1)
                        end
                %}
                
            end %field boundary loop
        end %direction loop

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %RESAVE POSITION TO ACCOUNT FOR CIRC SHIFT....
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    positionPerTimeStep=positionPerTime;
hasLapFieldInterference;

    %figure; histogram(abs(fieldWidths(2,:)),0:0.01:1)
    
    %try
save(unitDataFilePath,'hasLapFieldInterference','positionPerTimeStep','medianPos','isCircMaze', 'unwrappedSpeedPerTimeStep', 'numLapsPerDir','numFieldsPerDir','lapStartTimesPerDir',...
    'lapStopTimesPerDir','entryTimeIdxesPerDirPerFieldPerLap','exitTimeIdxesPerDirPerFieldPerLap','entryTimePerDirPerFieldPerLap',...
    'exitTimePerDirPerFieldPerLap','entryPosPerDirPerFieldPerLap','exitPosPerDirPerFieldPerLap','filledSignedSpeedPerTime','-append')
updatedData=load(unitDataFilePath)
    %catch ME
    %    disp(ME.message)
    %    disp('skipping...')
    %end
        
        

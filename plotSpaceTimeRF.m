function [] = plotSpaceTimeRF(spaceBins,timeBins,avgPhasePerSpaceTimeBin,fH)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    surfMode=0;
    polarMode=1;
    
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
        numPhaseBins=20;
        numPhaseBins=10;
        %numPhaseBins=8;
        phaseEdges=linspace(0,360,numPhaseBins+1);
        phaseBinAssignments=discretize(avgPhasePerSpaceTimeBin,phaseEdges);
        %phaseAndSpaceToTime=NaN(length(spaceBins),numPhaseBins);
        timesAccumulationHappens=0;
        phaseAndSpaceToTime=NaN(length(timeBins),numPhaseBins);
        phaseAndSpaceToTimeCount=zeros(length(timeBins),numPhaseBins);
        for phaseI=1:numPhaseBins
            spaceTimesPerPhase=find(phaseBinAssignments==phaseI);
            
            %associatedRows=NaN(length(spaceTimesPerPhase),1);
            %associatedCols=NaN(length(spaceTimesPerPhase),1);
            for i=1:length(spaceTimesPerPhase)
                [row,col]=ind2sub(size(avgPhasePerSpaceTimeBin),spaceTimesPerPhase(i));
                %if(isnan(phaseAndSpaceToTime(row,phaseI)))
                if(isnan(phaseAndSpaceToTime(col,phaseI)))
                    %phaseAndSpaceToTime(row,phaseI)=col;
                    phaseAndSpaceToTime(col,phaseI)=row;
                else
                    %phaseAndSpaceToTime(row,phaseI)=nanmean([phaseAndSpaceToTime(row,phaseI); col]);
                    phaseAndSpaceToTime(col,phaseI)=phaseAndSpaceToTime(col,phaseI)+row;
                    phaseAndSpaceToTimeCount(col,phaseI)=phaseAndSpaceToTimeCount(col,phaseI)+1;

                    timesAccumulationHappens=timesAccumulationHappens+1;
                end
                %associatedRows(i)=row;
                %associatedCols(i)=col;
            end
            %meanAssociatedRow=nanmean(associatedRows);
            %meanAssociatedCol=nanmean(associatedCols);
            %phaseAndSpaceToTime(round(meanAssociatedCol),phaseI)=meanAssociatedRow;
            
        end
        phaseAndSpaceToTimeCount(phaseAndSpaceToTimeCount==0)=NaN;
        phaseAndSpaceToTime=phaseAndSpaceToTime./phaseAndSpaceToTimeCount;
        polarplot3d(phaseAndSpaceToTime,'PlotType','surf')
        colormap(copper)
        timesAccumulationHappens
    else
        omarPcolor(spaceBins,timeBins,avgPhasePerSpaceTimeBin',fH)
        colormap(hsv)
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
        ylabel(cb,'Position')
    else
        ylabel(cb,'First spike theta phase (deg)')
    end
    %xlim([0 1])
    %ylim([-0.5 max(timeBins)+0.5])
end


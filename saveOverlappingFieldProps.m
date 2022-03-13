close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=1;
setTightSubplots_SpaceTime
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014', 'Cicero_09172014','Gatsby_08282013'};

%Cicero_09172014 has only 1 strongly precessing cell used
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

%halfLineHeight=0.01; %raster line height

    overlapFracThresh=0.4;

for si=1:length(sessionNames)
    %try
    currSesName=sessionNames{si};
    trackFilePath=sprintf('%sTrackProps.mat',currSesName);
    currTrackData=load(trackFilePath);
    positionTimeAxis=currTrackData.positionTimeAxis;
    
    rowColInSeqMatrix=currTrackData.rowColInSequenceMatrix;
    
    %currTrackData=currTrackData.trackPopInfo
    


    for di=1:2
        if(di==1)
           currDirStr='rightward';
        else
           currDirStr='leftward';
        end
        if(~isfield(rowColInSeqMatrix,currDirStr))
            continue
        end
        currOverlapMatrix=rowColInSeqMatrix.(currDirStr);
        unitNames=fieldnames(currTrackData.trackPopInfo.(currDirStr)); %RELYING ON STABLE NUMBERING OF FIELDNAMES
        
        for fi=1:size(currOverlapMatrix,1)
            for fii=1:size(currOverlapMatrix,2)
                
                if(currOverlapMatrix(fi,fii)==0)
                    continue
                end
                
                currUnitName=unitNames{fi};
                 overlappingNextUnitName=unitNames{fii};
                 
                currUnitNum=str2num(getSubstrBtwnSubstrs(currUnitName,'unit','Fi'));
                currFieldNum=str2num(currUnitName(end));
                fileBaseName=sprintf('%s_6-12Theta_Unit%dInfo',currSesName,currUnitNum);
                currUnitDataFilePath=fullfile(dataDir,sprintf('%s.mat',fileBaseName));
                currUnitData=load(currUnitDataFilePath);
                 if(~isfield(currUnitData.directionSpecificStats,currDirStr))
                    continue
                 end
                 
                currUnitStruct=currTrackData.trackPopInfo.(currDirStr).(currUnitName);
                nextUnitStruct=currTrackData.trackPopInfo.(currDirStr).(overlappingNextUnitName);
                
                [percentOverlap,similarFieldWidths]=getPercentOverlap(currUnitStruct,nextUnitStruct);
                
                if(percentOverlap<overlapFracThresh || ~similarFieldWidths)
                    continue
                end
                    
                    
                
                currUnitPhasePerTime=currUnitData.directionSpecificStats.(currDirStr).localSpikePhasePerTime;
                currUnitInFieldTimeIdxes=(currFieldNum==(currUnitData.directionSpecificStats.(currDirStr).fieldNumPerTime));
                
                localTimeInFieldPerTime=currUnitData.directionSpecificStats.(currDirStr).localTimeInFieldPerTime;
               
                
               
                
                nextUnitNum=str2num(getSubstrBtwnSubstrs(overlappingNextUnitName,'unit','Fi'));
                nextFieldNum=str2num(overlappingNextUnitName(end));
                nextUnitDataFilePath=fullfile(dataDir,sprintf('%s_6-12Theta_Unit%dInfo.mat',currSesName,nextUnitNum));
                nextUnitData=load(nextUnitDataFilePath);
                if(~isfield(nextUnitData.directionSpecificStats,currDirStr))
                    continue
                end
                nextUnitPhasePerTime=nextUnitData.directionSpecificStats.(currDirStr).localSpikePhasePerTime;
                nextUnitInFieldTimeIdxes=(nextFieldNum==(nextUnitData.directionSpecificStats.(currDirStr).fieldNumPerTime));
                
                
                inBothFieldsTimeIdxes=currUnitInFieldTimeIdxes & nextUnitInFieldTimeIdxes;
                currUnitPhasePerTime(~inBothFieldsTimeIdxes)=NaN;
                nextUnitPhasePerTime(~inBothFieldsTimeIdxes)=NaN;
                localTimeInFieldPerTime(~inBothFieldsTimeIdxes)=NaN;
                
                %phaseDiffPerTime=angdiffDeg([currUnitPhasePerTime(:)'; nextUnitPhasePerTime(:)']);
                  phaseDiffPerTime=angdiffDeg([nextUnitPhasePerTime(:)'; currUnitPhasePerTime(:)']);
                diffPhaseDiffPerTime=[NaN;diff(phaseDiffPerTime(:))];
                phaseDiffPerTime(diffPhaseDiffPerTime==0)=NaN; %just 1 per theta cycle
                figure(101); 
                subplot(3,1,1)
                plot(localTimeInFieldPerTime,phaseDiffPerTime,'k.','MarkerSize',20)
                hold on
                xlim([0 1])
                %ylim([0 180])
                ylim([-180 0])
                xlabel('Time in 1st field (frac)')
                ylabel('phase difference (deg)')
                daspect([1 360 1])
                
                
                title({removeUnderscores(sprintf('%s phase difference %s Vs %s',currSesName,currUnitName,overlappingNextUnitName)),sprintf('Percent overlap: %.2f',percentOverlap)})
              
                
                subplot(3,1,2)
                showFieldImage(currUnitData,di)
                title(sprintf('first field, %.2fm to %.2fm',currUnitStruct.fieldPosStart,currUnitStruct.fieldPosEnd))
                
                 subplot(3,1,3)
                 showFieldImage(nextUnitData,di)
                   title({sprintf('next field, %.2fm to %.2fm',nextUnitStruct.fieldPosStart,nextUnitStruct.fieldPosEnd)})

                   setFigFontTo(16)
                   maxFigMukkaalWidth
                saveas(gcf,sprintf('%s_PhaseDifference_%sVs%s.tif',currSesName,currUnitName,overlappingNextUnitName))
                
                close all
                
              
                %saveFileBaseName=sprintf('spacetimeCrossCorrelogram%sDir%d_%sVs%s',currSesName,di,currUnitName,overlappingNextUnitName);
                
            end
        end
    end
end

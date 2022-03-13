close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
showPlots=0;
setTightSubplots_SpaceTime
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014', 'Cicero_09172014','Gatsby_08282013'};

%Cicero_09172014 has only 1 strongly precessing cell used
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

%sessionNames={'Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Gatsby_08282013'};

minCenterSeparationM=0.1; %counting similar place fields is like double counting, dilutes effect
%minCenterSeparationM=0.2; %counting similar place fields is like double counting, dilutes effect
minFracOverlap=0.2;
minFracOverlap=0;

useEdgeConstraints=0;

for si=1:length(sessionNames)
    %try
    currSesName=sessionNames{si};
    
    [refUnitDataLeftward,refUnitDataRightward] = getRefUnitData(currSesName);
    
    trackFilePath=sprintf('%sTrackProps.mat',currSesName);
    currTrackData=load(trackFilePath);
    positionTimeAxis=currTrackData.positionTimeAxis;
    
    rowColInSeqMatrix=currTrackData.rowColInSequenceMatrix;
    minNumSpikesXcorr=30;
    
    
    %currTrackData=currTrackData.trackPopInfo

    for di=1:2
        if(di==1)
           currDirStr='rightward';
           refUnitData=refUnitDataRightward;
        else
           currDirStr='leftward';
           refUnitData=refUnitDataLeftward;
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
                
                currUnitStruct=currTrackData.trackPopInfo.(currDirStr).(currUnitName);
                nextUnitStruct=currTrackData.trackPopInfo.(currDirStr).(overlappingNextUnitName);
                
                [percentOverlap,similarFieldWidths]=getPercentOverlap(currUnitStruct,nextUnitStruct);
                
                if(percentOverlap<minFracOverlap || ~similarFieldWidths)
                    continue
                end
                currPairCenterSeparationM=abs(currUnitStruct.fieldPosCenterM - nextUnitStruct.fieldPosCenterM);
                currPairStartSeparationM=abs(currUnitStruct.fieldPosStart - nextUnitStruct.fieldPosStart);
                currPairEndSeparationM=abs(currUnitStruct.fieldPosEnd - nextUnitStruct.fieldPosEnd);
                
                %numLaps=refUnitData.directionSpecificStats.(currDirStr).
                %avgSpeedInMetaFieldPerLap=NaN(
                
                if(currPairCenterSeparationM<minCenterSeparationM)
                    continue
                end
                
                
                if(useEdgeConstraints)
                 if(currPairStartSeparationM<minCenterSeparationM || currPairEndSeparationM < minCenterSeparationM)
                    continue
                 end
                end
                
                saveFileBaseName=sprintf('TimeCrossCorrelogram%sDir%d_%sVs%s',currSesName,di,currUnitName,overlappingNextUnitName);
                
                if(exist([saveFileBaseName '.png'],'file'))
                    %continue
                end
                %try
                    fH=figure
                    subplot(2,2,[1 3])
                    %[currXcorrInfo]=getXcorrInfo(currUnitStruct,nextUnitStruct,positionTimeAxis,fH);
                    [currXcorrInfo]=getXcorrTimeInfo(currUnitStruct,nextUnitStruct,positionTimeAxis,refUnitData,fH);
                    
                    unit1OriginalNumSpikes=length(currUnitStruct.inFieldSpikeTimes);
                    unit2OriginalNumSpikes=length(nextUnitStruct.inFieldSpikeTimes);
                    
                    unit1NumSpikesPerSpeed=currXcorrInfo.numSpikesField1PerSpeed;
                    unit2NumSpikesPerSpeed=currXcorrInfo.numSpikesField2PerSpeed;
                    
                    if(min([unit1NumSpikesPerSpeed(:); unit2NumSpikesPerSpeed(:)])<minNumSpikesXcorr)
                        close all
                        continue
                    end
                    
                     title({removeUnderscores(sprintf('%s, %s, %s vs %s',currSesName,currDirStr,currUnitName,overlappingNextUnitName)),sprintf('unit 1 has %d (%d) low (high) speed spikes out of %d, unit 2 has %d (%d) low (high) speed spikes out of %d',...
                         unit1NumSpikesPerSpeed(1),unit1NumSpikesPerSpeed(2),unit1OriginalNumSpikes,unit2NumSpikesPerSpeed(1),unit2NumSpikesPerSpeed(2),unit2OriginalNumSpikes)});
                     axis square
                     subplot(2,2,2)  
                      imshow(currUnitStruct.imagePath)
                      subplot(2,2,4)  
                       imshow(nextUnitStruct.imagePath)

                        setFigFontTo(18)
                         maxFig
                         print('-r75',saveFileBaseName,'-dpng')
                        close all

                        eval(sprintf('%sVs%sTimeXcorrInfoWithBehavTimes = currXcorrInfo;',currUnitName,overlappingNextUnitName))

                    save(trackFilePath,sprintf('%sVs%sTimeXcorrInfoWithBehavTimes',currUnitName,overlappingNextUnitName),'-append')
                  %{  
                catch ME
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    disp('SKIPPING:')
                    disp(ME.message)
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

                end
                    %}
            end
        end
    end
    %{
      catch ME
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                    disp('SKIPPING:')
                    disp(ME.message)
                    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

                end
    %}
end
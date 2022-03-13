
close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';


%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014', 'Cicero_09172014','Gatsby_08282013'};

%Cicero_09172014 has only 1 strongly precessing cell used
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013', 'Cicero_09102014','Cicero_09172014','Gatsby_08282013'};
%sessionNames={'Achilles_11012013'};

minSpeed=0.05; %m/s

lapFig=figure(10101)

for si=1:length(sessionNames)
    currSesName=sessionNames{si};
    currSesFilePaths=getFilePathsRegex(dataDir,sprintf('%s*mat',currSesName));
    
    [leftwardRefUnitData,rightwardRefUnitData]=getRefUnitData(currSesName);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get spike raster of fields sorted by position center
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trackPopInfo=[];
    for fi=1:length(currSesFilePaths)
        currFilePath=currSesFilePaths{fi};
        currUnitData=load(currFilePath);
        currUnitID=currUnitData.unitInfo.unitIDnum;
        imagePaths=currUnitData.imagePaths;
        
        
        
        positionTimeAxis=currUnitData.positionTimeAxis;
        filledSignedSpeedPerTime=currUnitData.filledSignedSpeedPerTime;
        
        for di=1:2
     
            if(di==1)
                dirStr='right';
                refUnitData=rightwardRefUnitData;
                if(~isfield(currUnitData,'manualFieldEndsMrightward'))
                    continue
                end
                manualFieldEndsM=currUnitData.manualFieldEndsMrightward;
                manualFieldStartsM=currUnitData.manualFieldStartsMrightward;
            else
                dirStr='left';
                refUnitData=leftwardRefUnitData;
                 if(~isfield(currUnitData,'manualFieldEndsMleftward'))
                    continue
                end
                manualFieldEndsM=currUnitData.manualFieldEndsMleftward;
                manualFieldStartsM=currUnitData.manualFieldStartsMleftward;
            end
            currDirStr=[dirStr 'ward'];
            
            if(~isfield(currUnitData.directionSpecificStats,currDirStr))
                continue
            end
                
            
            if(~iscell(imagePaths))
                imagePath=imagePaths;
            else
                if(contains(imagePaths{1},currDirStr))
                    imagePath=imagePaths{1};
                else
                    imagePath=imagePaths{2};
                end
            end
            
            currNumFields=length(manualFieldStartsM);
            
            for fii=1:currNumFields
                currFieldIDstr=sprintf('unit%dField%d',currUnitID,fii);
                currFieldStartPos=manualFieldStartsM(fii);
                currFieldEndPos=manualFieldEndsM(fii);
                
                if(isnan(manualFieldStartsM(fii)) || isnan(manualFieldEndsM(fii)))
                    continue
                end
                
                currDirSpikePositions=currUnitData.(sprintf('%sSpikePositions',dirStr));
                currDirSpikeTimes=currUnitData.(sprintf('%sSpikeTimes',dirStr));
                currDirSpikePhases=currUnitData.(sprintf('%sSpikePhases',dirStr));
                
                lapNumPerTime=currUnitData.directionSpecificStats.(currDirStr).lapNumberPerTime;
                
                
                %lapNumPerTime=refUnitData.directionSpecificStats.(currDirStr).lapNumberPerTime;

                
                figure(lapFig)
                allFieldsMinTime=19648-290;
                plot(currUnitData.positionTimeAxis-allFieldsMinTime,lapNumPerTime,'LineWidth',5)
                hold on
                currDirSpikeLapNums=interp1(positionTimeAxis,lapNumPerTime,currDirSpikeTimes);
                
                currDirSpikeSpeeds=interp1(positionTimeAxis,filledSignedSpeedPerTime,currDirSpikeTimes);
         
                
                %for leftward direction, end position < start position
                if(di==1)
                    inFieldSpikeIdxes=currDirSpikePositions<=currFieldEndPos & currDirSpikePositions>=currFieldStartPos;
                      goodSpeedSpikeIdxes=currDirSpikeSpeeds>=minSpeed; 
                else
                    inFieldSpikeIdxes=currDirSpikePositions>=currFieldEndPos & currDirSpikePositions<=currFieldStartPos;   
                    goodSpeedSpikeIdxes=currDirSpikeSpeeds<=(-minSpeed); 
                end
              
                   
                inFieldSpikeTimes=currDirSpikeTimes(inFieldSpikeIdxes);  
                
                inFieldSpikeLapNums=currDirSpikeLapNums(inFieldSpikeIdxes);
                inFieldSpikeLapNumsGoodSpeed=currDirSpikeLapNums(inFieldSpikeIdxes & goodSpeedSpikeIdxes); 
                
                inFieldSpikeTimesGoodSpeed=currDirSpikeTimes(inFieldSpikeIdxes & goodSpeedSpikeIdxes); 
                
                
                inFieldSpikePhases=currDirSpikePhases(inFieldSpikeIdxes);  
                inFieldSpikePhasesGoodSpeed=currDirSpikePhases(inFieldSpikeIdxes & goodSpeedSpikeIdxes); 
                
                inFieldSpikePositions=currDirSpikePositions(inFieldSpikeIdxes);  
                inFieldSpikePositionsGoodSpeed=currDirSpikePositions(inFieldSpikeIdxes & goodSpeedSpikeIdxes); 
                
                inFieldSpikeSpeeds=currDirSpikeSpeeds(inFieldSpikeIdxes);
                
                fieldPosCenterM=mean([currFieldStartPos currFieldEndPos]);
                
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikeTimes=inFieldSpikeTimes;
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikePhases=inFieldSpikePhases;
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikePositions=inFieldSpikePositions;
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikeSpeeds=inFieldSpikeSpeeds;
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikeLapNums=inFieldSpikeLapNums;
                
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikeTimesGoodSpeed=inFieldSpikeTimesGoodSpeed;
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikePhasesGoodSpeed=inFieldSpikePhasesGoodSpeed;
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikePositionsGoodSpeed=inFieldSpikePositionsGoodSpeed;
                trackPopInfo.(currDirStr).(currFieldIDstr).inFieldSpikeLapNumsGoodSpeed=inFieldSpikeLapNumsGoodSpeed;
                
                trackPopInfo.(currDirStr).(currFieldIDstr).fieldPosCenterM=fieldPosCenterM;  
                trackPopInfo.(currDirStr).(currFieldIDstr).fieldPosStart=currFieldStartPos;  
                trackPopInfo.(currDirStr).(currFieldIDstr).fieldPosEnd=currFieldEndPos;  
                trackPopInfo.(currDirStr).(currFieldIDstr).imagePath=imagePath;  
            end
        end
        %trackPopInfo
    end
    si
    save(sprintf('%sTrackProps.mat',currSesName),'trackPopInfo','positionTimeAxis')
    
end
close all; clear all; clc

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
timeBoundSaveDir='./hc3Images/spaceTimePhase';

touchDir(timeBoundSaveDir)

unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

unitFilePaths=getFilePathsRegex(unitDataDir,'*mat');

totalFieldCount=0;
%setTightSubplots

%redoBounds=1;

minRunSpeedDisp=0.05;%m/s
maxRunSpeedDisp=0.8; %m/s

plotByLap=0; 
plotBySpeed=0;
plotByDist=1;

startFi=1;
%startFi=62;
startFi=200
for fi=startFi:length(unitFilePaths)
    resaveImgFlag=0;
    
    if(mod(fi,50)==0)
        disp(fi)
    end
    currUnitFilePath=unitFilePaths{fi};
    currUnitStruct=load(currUnitFilePath);

        currUnitData=currUnitStruct.unitInfo;

        currUnitFileName=getFileNameFromPath(currUnitFilePath);
        
        cycleTimeFieldBoundPerDirPerField=currUnitStruct.cycleTimeFieldBoundPerDirPerField;

        numRightwardFields=length(currUnitStruct.rightwardFieldStartEndM)/2;
        numLeftwardFields=length(currUnitStruct.leftwardFieldStartEndM)/2;

        maxNumFieldsPerDir=max([numRightwardFields numLeftwardFields]);

        %if(redoBounds)
            rightwardTimeBounds=NaN(maxNumFieldsPerDir,1);
        %end
        
        for ri=1:numRightwardFields
            currFieldThetaData=currUnitStruct.thetaBasedTimeVsPhaseInfo.('rightward').(sprintf('field%d',ri));
            baseName=sprintf('%s_Cell%d_Rightward_Field%d',currUnitData.sessionName,currUnitData.unitIDnum,ri);
            
            saveImgFilePath=fullfile(timeBoundSaveDir,sprintf('%s_Cell%dRightwardField%d_SpaceTimePhase.tif',currUnitData.sessionName,currUnitData.unitIDnum,ri));
            
            spikeLapNums=currFieldThetaData.allLapInFieldSpikeLapNums;
            spikeSpeeds=abs(currFieldThetaData.allLapInFieldSpikeSpeeds);
            spikeDistsFrac=currFieldThetaData.allLapInFieldSpikeDistsFieldFrac;
            
            currTimeBound=cycleTimeFieldBoundPerDirPerField{ri,1};
            
          
            
            %lapColors=jet(max(spikeLapNums));

            fH=figure; 
          
            spikeTimeFieldFrac=currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry/currTimeBound;
          scatter(spikeDistsFrac,spikeTimeFieldFrac,50, currFieldThetaData.allLapInFieldSpikePhases,'filled')
            colormap(jet)
            cb=colorbar
            hold on
            ylabel(cb,'Theta phase')
            caxis([0 360])
    
            xlim([0 1])
            ylim([0 1])
               
            ylabel('Time in field (frac)')
            xlabel('Dist in field (frac)')
            setFigFontTo(18)
            
            hold on
            
            set(gca,'Color','k')
            
            
         
                set(gcf, 'InvertHardcopy', 'off')
                %rightwardTimeBounds(ri)=currTimeBound;
               
                %maxFig
                 axis square
                 uberTitle(removeUnderscores(baseName))
                saveas(gcf,saveImgFilePath)
                close(fH)
            
            
        end

        leftwardTimeBounds=NaN(maxNumFieldsPerDir,1);
        for ri=1:numLeftwardFields
            currFieldThetaData=currUnitStruct.thetaBasedTimeVsPhaseInfo.('leftward').(sprintf('field%d',ri));
            saveImgFilePath=fullfile(timeBoundSaveDir,sprintf('%s_Cell%dLeftwardField%d_SpaceTimePhase.tif',currUnitData.sessionName,currUnitData.unitIDnum,ri));
            baseName=sprintf('%s_Cell%d_Leftward_Field%d',currUnitData.sessionName,currUnitData.unitIDnum,ri);
            spikeLapNums=currFieldThetaData.allLapInFieldSpikeLapNums;
            spikeSpeeds=abs(currFieldThetaData.allLapInFieldSpikeSpeeds);
            spikeDistsFrac=currFieldThetaData.allLapInFieldSpikeDistsFieldFrac;
            currTimeBound=cycleTimeFieldBoundPerDirPerField{ri,2};
            
            
            
            %lapColors=jet(max(spikeLapNums));

            fH=figure; 
            spikeTimeFieldFrac=currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry/currTimeBound;
          scatter(spikeDistsFrac,spikeTimeFieldFrac,50, currFieldThetaData.allLapInFieldSpikePhases,'filled')
            colormap(jet)
            cb=colorbar
            hold on
            ylabel(cb,'Theta phase')
            caxis([0 360])
    
            xlim([0 1])
            ylim([0 1])
               
            ylabel('Time in field (frac)')
            xlabel('Dist in field (frac)')
            
            hold on
            
            set(gca,'Color','k')
            
           
            
         
                set(gcf, 'InvertHardcopy', 'off')
                leftwardTimeBounds(ri)=currTimeBound;
                axis square
                %maxFig
                 uberTitle(removeUnderscores(baseName))
                saveas(gcf,saveImgFilePath)
                close(fH)
            
            
         
        end

    
end
%totalFieldCount/2

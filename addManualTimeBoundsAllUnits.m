close all; clear all; clc

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
timeBoundSaveDir='./hc3Images/timeVsPhaseWithBound';

touchDir(timeBoundSaveDir)

unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

unitFilePaths=getFilePathsRegex(unitDataDir,'*mat');

totalFieldCount=0;
%setTightSubplots
redoBounds=1;
redoBounds=0;
%redoBounds=1;

minRunSpeedDisp=0.05;%m/s
maxRunSpeedDisp=0.8; %m/s

plotByLap=0; 
plotBySpeed=0;
plotByDist=1;

for fi=1:length(unitFilePaths)
    resaveImgFlag=0;
    
    if(mod(fi,50)==0)
        disp(fi)
    end
    currUnitFilePath=unitFilePaths{fi};
    currUnitStruct=load(currUnitFilePath);
    
    if(isfield(currUnitStruct,'cycleTimeFieldBoundPerDirPerField'))
        currUnitData=currUnitStruct.unitInfo;

        currUnitFileName=getFileNameFromPath(currUnitFilePath);
        
        cycleTimeFieldBoundPerDirPerField=currUnitStruct.cycleTimeFieldBoundPerDirPerField;

        numRightwardFields=length(currUnitStruct.rightwardFieldStartEndM)/2;
        numLeftwardFields=length(currUnitStruct.leftwardFieldStartEndM)/2;

        maxNumFieldsPerDir=max([numRightwardFields numLeftwardFields]);

        rightwardTimeBounds=NaN(maxNumFieldsPerDir,1);
        for ri=1:numRightwardFields
            currFieldThetaData=currUnitStruct.thetaBasedTimeVsPhaseInfo.('rightward').(sprintf('field%d',ri));
            baseName=sprintf('%s_Cell%d_Rightward_Field%d',currUnitData.sessionName,currUnitData.unitIDnum,ri);
            
            saveImgFilePath=fullfile(timeBoundSaveDir,sprintf('%s_Cell%dRightwardField%d_CycleTimeVsPhase.tif',currUnitData.sessionName,currUnitData.unitIDnum,ri));
            
            spikeLapNums=currFieldThetaData.allLapInFieldSpikeLapNums;
            spikeSpeeds=abs(currFieldThetaData.allLapInFieldSpikeSpeeds);
            spikeDistsFrac=currFieldThetaData.allLapInFieldSpikeDistsFieldFrac;
            
            currTimeBound=cycleTimeFieldBoundPerDirPerField{ri,1};
            
            if(redoBounds)
                if(~isnan(currTimeBound))
                    rightwardTimeBounds(ri)=currTimeBound;
                    continue
                end
            end
            
            %lapColors=jet(max(spikeLapNums));

            fH=figure; 
            if(plotByLap)
                scatter(currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry,currFieldThetaData.allLapInFieldSpikePhases,50,spikeLapNums,'filled')
            elseif(plotBySpeed)
                scatter(currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry,currFieldThetaData.allLapInFieldSpikePhases,50,spikeSpeeds,'filled')
            elseif(plotByDist)
                scatter(currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry,currFieldThetaData.allLapInFieldSpikePhases,50,spikeDistsFrac,'filled')
            end
            colormap(jet)
            cb=colorbar
            hold on
         
            
           if(plotByLap)
                ylabel(cb,'Lap number')
            elseif(plotBySpeed)
                ylabel(cb,'Running speed (m/s)')
                caxis([minRunSpeedDisp maxRunSpeedDisp])
            elseif(plotByDist)
                ylabel(cb,'Dist in field (fraction)')
                caxis([0 1])
            end
            xlim([0 30])
            ylim([0 360])
            hold on
            
               plot([currTimeBound currTimeBound],ylim,'w--','LineWidth',5)
               
            xlabel('No. of theta cycles in field')
            ylabel('Theta phase')
            setFigFontTo(18)
            
            hold on
            
            set(gca,'Color','k')
            uberTitle(removeUnderscores(baseName))
            
            fi
            if(redoBounds)
                [currFieldTimeBoundNew,~] = ginput

                if(length(currFieldTimeBoundNew)==1)
                    rightwardTimeBounds(ri)=currFieldTimeBoundNew;
                else
                    rightwardTimeBounds(ri)=currTimeBound;
                end
                close(fH)
            else
                set(gcf, 'InvertHardcopy', 'off')
                rightwardTimeBounds(ri)=currTimeBound;
                saveas(gcf,saveImgFilePath)
                close(fH)
            end
            
        end

        leftwardTimeBounds=NaN(maxNumFieldsPerDir,1);
        for ri=1:numLeftwardFields
            currFieldThetaData=currUnitStruct.thetaBasedTimeVsPhaseInfo.('leftward').(sprintf('field%d',ri));
            saveImgFilePath=fullfile(timeBoundSaveDir,sprintf('%s_Cell%dLeftwardField%d_CycleTimeVsPhase.tif',currUnitData.sessionName,currUnitData.unitIDnum,ri));
            baseName=sprintf('%s_Cell%d_Leftward_Field%d',currUnitData.sessionName,currUnitData.unitIDnum,ri);
            spikeLapNums=currFieldThetaData.allLapInFieldSpikeLapNums;
            spikeSpeeds=abs(currFieldThetaData.allLapInFieldSpikeSpeeds);
            spikeDistsFrac=currFieldThetaData.allLapInFieldSpikeDistsFieldFrac;
            currTimeBound=cycleTimeFieldBoundPerDirPerField{ri,2};
            
            if(redoBounds)
             if(~isnan(currTimeBound))
                leftwardTimeBounds(ri)=currTimeBound;
                continue
             end
            end
            
            
            %lapColors=jet(max(spikeLapNums));

            fH=figure; 
            if(plotByLap)
                scatter(currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry,currFieldThetaData.allLapInFieldSpikePhases,50,spikeLapNums,'filled')
            elseif(plotBySpeed)
                scatter(currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry,currFieldThetaData.allLapInFieldSpikePhases,50,spikeSpeeds,'filled')
            elseif(plotByDist)
                scatter(currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry,currFieldThetaData.allLapInFieldSpikePhases,50,spikeDistsFrac,'filled')
            end
            colormap(jet)
            cb=colorbar
            hold on
         
            if(plotByLap)
                ylabel(cb,'Lap number')
            elseif(plotBySpeed)
                ylabel(cb,'Running speed (m/s)')
                caxis([minRunSpeedDisp maxRunSpeedDisp])
            elseif(plotByDist)
                ylabel(cb,'Dist in field (fraction)')
                caxis([0 1])
            end
            xlim([0 30])
            ylim([0 360])
            
               plot([currTimeBound currTimeBound],ylim,'w--','LineWidth',5)
               
            xlabel('No. of theta cycles in field')
            ylabel('Theta phase')
            setFigFontTo(18)
            
            hold on
            
            set(gca,'Color','k')
            
            uberTitle(removeUnderscores(baseName))
            
            fi
            if(redoBounds)
                [currFieldTimeBoundNew,~] = ginput

                if(length(currFieldTimeBoundNew)==1)
                    leftwardTimeBounds(ri)=currFieldTimeBoundNew;
                else
                    leftwardTimeBounds(ri)=currTimeBound;
                end
                close(fH)
            else
                set(gcf, 'InvertHardcopy', 'off')
                leftwardTimeBounds(ri)=currTimeBound;
                saveas(gcf,saveImgFilePath)
                close(fH)
            end
            
         
        end
    else
        currUnitData=currUnitStruct.unitInfo;

        currUnitFileName=getFileNameFromPath(currUnitFilePath);

        numRightwardFields=length(currUnitStruct.rightwardFieldStartEndM)/2;
        numLeftwardFields=length(currUnitStruct.leftwardFieldStartEndM)/2;

        maxNumFieldsPerDir=max([numRightwardFields numLeftwardFields]);

        rightwardTimeBounds=NaN(maxNumFieldsPerDir,1);
        for ri=1:numRightwardFields
            currFieldThetaData=currUnitStruct.thetaBasedTimeVsPhaseInfo.('rightward').(sprintf('field%d',ri));

            fH=figure; 
            plot(currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry,currFieldThetaData.allLapInFieldSpikePhases,'k.','MarkerSize',15)
            xlim([0 30])
            ylim([0 360])
            xlabel('No. of theta cycles in field')
            ylabel('Theta phase')
            setFigFontTo(18)


            [currFieldTimeBound,~] = ginput

            %assert(length(currFieldTimeBound)==1)

             if(length(currFieldTimeBound)==1)
                rightwardTimeBounds(ri)=currFieldTimeBound;
             end
            close(fH)
        end

        leftwardTimeBounds=NaN(maxNumFieldsPerDir,1);
        for ri=1:numLeftwardFields
            currFieldThetaData=currUnitStruct.thetaBasedTimeVsPhaseInfo.('leftward').(sprintf('field%d',ri));

            fH=figure; 
            plot(currFieldThetaData.allLapSpikeCycleIDsFromFieldEntry,currFieldThetaData.allLapInFieldSpikePhases,'k.','MarkerSize',15)
            xlim([0 30])
            ylim([0 360])
            xlabel('No. of theta cycles in field')
            ylabel('Theta phase')
            setFigFontTo(18)

            [currFieldTimeBound,~] = ginput

            if(length(currFieldTimeBound)==1)
                leftwardTimeBounds(ri)=currFieldTimeBound;
            end
            close(fH)
        end


       
    end
    
     cycleTimeFieldBoundPerDirPerField=table(rightwardTimeBounds,leftwardTimeBounds,'VariableNames',{'rightward','leftward'})

        save(currUnitFilePath,'cycleTimeFieldBoundPerDirPerField', '-append')
    
end
%totalFieldCount/2

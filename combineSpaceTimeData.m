close all
dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir';

trackDir=1;

filePaths=getFilePathsRegex(dataDir,sprintf('*dir%d*mat',trackDir));

numSpaceBins=50;
%numSpaceBins=25;
%numSpaceBins=15;
numTimeBins=15;
%numTimeBins=10;

trackLength=120; %cm

numFilesToUse=length(filePaths);
%numFilesToUse=100;
totalFieldResponseAcrossCells=NaN(numSpaceBins,numTimeBins,numFilesToUse);


%for i=1:length(filePaths)
%for i=1:9
skipCount=0;
for i=1:numFilesToUse
    
    currFilePath=filePaths{i};
    currData=load(currFilePath);
    
    spaceVals=currData.spaceTimePhaseInfo.spaceBins;
    if(isempty(spaceVals) || sum(isnan(spaceVals)) >0)
        continue
    end
    phaseVsSpace=currData.spaceTimePhaseInfo.avgPhasePerSpaceBin;
    phaseVsSpaceCov=corrcoef(spaceVals,phaseVsSpace);
    
    
    phaseVsSpaceR=phaseVsSpaceCov(1,2);
     if(abs(phaseVsSpaceR)<0)
        skipCount=skipCount+1;
        continue  
    end
    
    absSpaceTimePhaseTriplets=currData.spaceTimePhaseInfo.spaceTimePhasePerCycle;
    minSpaceFrac=min(absSpaceTimePhaseTriplets(:,1))/max(absSpaceTimePhaseTriplets(:,1));
    startSpaceFracTrack=min(currData.spikePerCycleInfo.peakFieldStartPosNorm(trackDir),minSpaceFrac );
    endSpaceFracTrack=currData.spikePerCycleInfo.peakFieldEndPosNorm(trackDir);
    fieldWidthFracTrack=endSpaceFracTrack-startSpaceFracTrack;
    
    spaceBinWidth=(endSpaceFracTrack-startSpaceFracTrack)/numSpaceBins;
    meanSpeed=currData.spikePerCycleInfo.meanSpeed;
    speedInLapsPerSecond=meanSpeed/trackLength;
   
    startTimeFracTrack=startSpaceFracTrack/speedInLapsPerSecond;
    endTimeFracTrack=endSpaceFracTrack/speedInLapsPerSecond;
    fieldDurationFracTrack=endTimeFracTrack-startTimeFracTrack;
    timeBinWidth=(endTimeFracTrack-startTimeFracTrack)/numTimeBins;
    
    fieldSpaceBinWidth=1/numSpaceBins;
    fieldTimeBinWidth=1/numTimeBins;
    
    %commonNormSpaceBinCenters=(spaceBinWidth/2):spaceBinWidth:(numSpaceBins-spaceBinWidth/2);
    fieldNormSpaceBinCenters=(fieldSpaceBinWidth/2):fieldSpaceBinWidth:(1-fieldSpaceBinWidth/2);
    %commonNormTimeBinCenters=(timeBinWidth/2):timeBinWidth:(numTimeBins-timeBinWidth/2);
    fieldNormTimeBinCenters=(fieldTimeBinWidth/2):fieldTimeBinWidth:(1-fieldTimeBinWidth/2);
    
 
    numDataCycles=size(absSpaceTimePhaseTriplets,1);
    
    commonFieldPhaseResponse=NaN(numSpaceBins,numTimeBins);
    
    spaceBinFact=40;
    spaceBinFact=max(absSpaceTimePhaseTriplets(:,1));
    %timeBinFact=meanSpeed; %mistake from getSpaceTimePhaseInfo, need to use again here
    timeBinFact=max(absSpaceTimePhaseTriplets(:,2));
    %figure
    for si=1:numSpaceBins
        for ti=1:numTimeBins
            isInCurrentBin=zeros(numDataCycles,1);
            wasHere=0;
            for ci=1:numDataCycles
                currPosFracTrack=absSpaceTimePhaseTriplets(ci,1)/spaceBinFact;
                currPosFracField=(currPosFracTrack-startSpaceFracTrack)/fieldWidthFracTrack;
                %si/numSpaceBins
                
                currTimeFracTrack=absSpaceTimePhaseTriplets(ci,2)/timeBinFact;
                %currTimeFracField=(currTimeFracTrack-startTimeFracTrack)/fieldDurationFracTrack;
                %currTimeFracField=(currTimeFracTrack)*fieldDurationFracTrack;
                currTimeFracField=(currTimeFracTrack);
                %ti/numTimeBins
                
                currMinPosFracField=(fieldNormSpaceBinCenters(si)-fieldSpaceBinWidth/2);
                currMaxPosFracField=(fieldNormSpaceBinCenters(si)+fieldSpaceBinWidth/2);
                currMinTimeFracField=(fieldNormTimeBinCenters(ti)-fieldTimeBinWidth/2);
                currMaxTimeFracField=(fieldNormTimeBinCenters(ti)+fieldTimeBinWidth/2);
                 if(false && mod(ci,10)==0) %%&& mod(ti,10)==0 && mod(si,10)==0)
                    plot([currMinPosFracField currMinPosFracField],ylim,'k-')
                    xlim([0 2])

                    hold on
                    plot([currMaxPosFracField currMaxPosFracField],ylim,'k-')
                    plot(xlim,[currMinTimeFracField currMinTimeFracField],'b-')
                    ylim([0 2])
                    plot(xlim,[currMaxTimeFracField currMaxTimeFracField],'b-')
                    plot(currPosFracField,currTimeFracField,'ro','MarkerSize',5)
                    title(sprintf('Curr pos: %.5f, Curr time: %.5f, Field start %.4f, Field end %.4f',currPosFracField,currTimeFracField, ...
                        startSpaceFracTrack,endTimeFracTrack))
              
                    drawnow
                    %pause(0.5)
                end
                if(currMinPosFracField<currPosFracField && currPosFracField<=currMaxPosFracField ...
                        && currMinTimeFracField<currTimeFracField && currTimeFracField<=currMaxTimeFracField)
                    isInCurrentBin(ci)=1;
                    %disp('here')
                    wasHere=1;
                end
            end
            if(wasHere==0)
                %fds
            end
            phasesInCurrentSpaceTimeBin=absSpaceTimePhaseTriplets(logical(isInCurrentBin),3);
            commonFieldPhaseResponse(si,ti)=circMeanDeg(phasesInCurrentSpaceTimeBin);
        end
        
    end
    totalFieldResponseAcrossCells(:,:,i)=commonFieldPhaseResponse;
    
    fH=figure
    subplot(2,1,1)
   
    plotSpaceTimeRF(currData.spaceTimePhaseInfo.spaceBins,currData.spaceTimePhaseInfo.timeBins,currData.spaceTimePhaseInfo.avgPhasePerSpaceTimeBin,fH)
    title('absolute spacetime')
    subplot(2,1,2)
    omarPcolor(fieldNormSpaceBinCenters,fieldNormTimeBinCenters,commonFieldPhaseResponse',fH)
    %omarPcolor(fieldNormSpaceBinCenters,fieldNormTimeBinCenters,totalFieldResponseAcrossCells',fH)

    title('Just field')
    currFileName=getFileNameFromPath(currFilePath);
    fileNameRoot=currFileName(1:(end-4));
    saveas(gcf,sprintf('%s.tif',fileNameRoot))
end
%%
figure
allCellMeanPhaseSpaceTime=circMeanDeg(totalFieldResponseAcrossCells,3);
omarPcolor(fieldNormSpaceBinCenters,fieldNormTimeBinCenters,allCellMeanPhaseSpaceTime')
colormap(hsv)
cb=colorbar

phaseCrossSectionsTime=NaN(numSpaceBins,numTimeBins);
%%
cmp=copper(numSpaceBins);
figure
for si=1:numSpaceBins
    plot(fieldNormTimeBinCenters,allCellMeanPhaseSpaceTime(si,:),'Color',cmp(si,:))
    hold on
    xlabel('Time')
end

figure
cmp=copper(numTimeBins);
for ti=1:numTimeBins
    plot(fieldNormSpaceBinCenters,allCellMeanPhaseSpaceTime(:,ti),'Color',cmp(ti,:))
    hold on
     xlabel('Distance within field')
end
skipCount
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load unit struct file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  clc

%NOT RELOADING LARGE FILE, TAKE OUT TO RELOAD
%clear all;
clearvars -except spikeDataPerField

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');

filePaths=getFilePathsRegex(unitDataDir,'*mat');

tic
disp('loading spike data....')
if(~exist('spikeDataPerField','var'))
    spikeDataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
end
toc

spikeTimesInExpPerField=spikeDataPerField.spikeTimesInExpPerField;
spikePhasesPerField =spikeDataPerField.spikePhasesPerField;

distInFieldPerSpikePerField=spikeDataPerField.distInFieldPerSpikePerField;
spikeTimeFracInFieldPerField=spikeDataPerField.spikeTimeFracInFieldPerField;
unitInfoPathPerField=spikeDataPerField.unitInfoPathPerField;

%spikeSpeedsPerField=spikeDataPerField.spikeSpeedsPerField;

spikeSpeedsPerField=spikeDataPerField.spikeAvgSpeedInTraversalsPerField;

accPerSpikePerField=spikeDataPerField.accPerSpikePerField;

minSpeed=0;
maxSpeed=Inf;
%maxAccMPerSecSquared=0.1; %10 cm/s^2
%maxAccMPerSecSquared=Inf; %10 cm/s^2
%minAccMPerSecSquared=-0.1; %10 cm/s^2
%minAccMPerSecSquared=-0.01; %10 cm/s^2
%minAccMPerSecSquared=-0.2; %10 cm/s^2 
%minAccMPerSecSquared=0; %10 cm/s^2 %4x

%%4x+
minAccMPerSecSquared=-0.1; %10 cm/s^2 
maxAccMPerSecSquared=0.1; %10 cm/s^2

%minAccMPerSecSquared=-0.05; %5 cm/s^2 
%maxAccMPerSecSquared=0.05; %5 cm/s^2

%minAccMPerSecSquared=-0.15; %15 cm/s^2 
%maxAccMPerSecSquared=0.15; %15 cm/s^2

%minAccMPerSecSquared=-Inf; %10 cm/s^2 
%maxAccMPerSecSquared=Inf; %10 cm/s^2
%minAccMPerSecSquared=-Inf; %3x
%maxSpeed=0.25;


totalFieldCount=spikeDataPerField.totalFieldCount;

showPlots=0;

startFieldIdx=1;
 fH=figure;
 
 behavioralBinEdges=linspace(0,1,10);


 phaseColors=jet(360);
 
 usedFieldCount=0;
 allDistTimePhaseTriples=[];
for fi=startFieldIdx:totalFieldCount
     tic

       fi
    currFilePath=unitInfoPathPerField{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));

    dataStruct=load(currFilePath);

    %currFieldSpikeTimesInExp=spikeTimesInExpPerField{fi};

    currFieldTimeInFieldFracPerSpike=spikeTimeFracInFieldPerField{fi};
    currFieldDistInFieldFracPerSpike=distInFieldPerSpikePerField{fi};
    currFieldSpikePhases=spikePhasesPerField{fi};
    
    currFieldSpikeSpeeds=spikeSpeedsPerField{fi};
    currFieldSpikeAccs=accPerSpikePerField{fi};
    
    
    if(~isempty(currFieldSpikeSpeeds))
         validSpeedIdxes=abs(currFieldSpikeSpeeds)>=minSpeed & abs(currFieldSpikeSpeeds)<=maxSpeed & currFieldSpikeAccs<=maxAccMPerSecSquared  & currFieldSpikeAccs>=minAccMPerSecSquared;
    else
         validSpeedIdxes=currFieldSpikeAccs<=maxAccMPerSecSquared & currFieldSpikeAccs>=minAccMPerSecSquared ;

    end
    
    currFieldTimeInFieldFracPerSpike=currFieldTimeInFieldFracPerSpike(validSpeedIdxes);
    currFieldDistInFieldFracPerSpike=currFieldDistInFieldFracPerSpike(validSpeedIdxes);
    currFieldSpikePhases=currFieldSpikePhases(validSpeedIdxes);
    
    if(isempty(currFieldSpikePhases))
        continue
    end
    
    
    [bestShift,bestShiftedPhases]=getBestShiftDeg(currFieldSpikePhases,currFieldDistInFieldFracPerSpike,currFieldTimeInFieldFracPerSpike);
    
     allDistTimePhaseTriples=[allDistTimePhaseTriples; [currFieldDistInFieldFracPerSpike(:), currFieldTimeInFieldFracPerSpike(:), bestShiftedPhases(:)]];

     usedFieldCount=usedFieldCount+1;
    %allDistTimePhaseTriples=[allDistTimePhaseTriples; [currFieldDistInFieldFracPerSpike(:), currFieldTimeInFieldFracPerSpike(:), currFieldSpikePhases(:)]];

    getBinAvgs=0;
    if(getBinAvgs)
        behavioralBinEdges
        [timeBinCenters,phaseAvgValuesPerTime,lowerLimPerTimeBin,upperLimPerTimeBin]=getBinnedCircAverages(currFieldTimeInFieldFracPerSpike,currFieldSpikePhases,behavioralBinEdges);
       [distBinCenters,phaseAvgValuesPerDist,lowerLimPerDistBin,upperLimPerDistBin]=getBinnedCircAverages(currFieldDistInFieldFracPerSpike,currFieldSpikePhases,behavioralBinEdges);
    end
    if(showPlots)
        figure(fH)
        %subplot(10,10,fi)
        
        %{
        for si=1:length(currFieldDistInFieldFracPerSpike)
            currSpikePhase=round(currFieldSpikePhases(si));
            
            if(currSpikePhase<1 || currSpikePhase>360)
                currSpikePhase=1;
            end

            plot(currFieldDistInFieldFracPerSpike(si),currFieldTimeInFieldFracPerSpike(si),'.','Color',phaseColors(currSpikePhase,:))
            hold on
        
        end
        %}
        
        
        %{
        errorbar(timeBinCenters,phaseAvgValuesPerTime,abs(angdiffDeg([upperLimPerTimeBin'; phaseAvgValuesPerTime'])),abs(angdiffDeg([lowerLimPerTimeBin'; phaseAvgValuesPerTime'])),'b')
        hold on
         errorbar(distBinCenters,phaseAvgValuesPerDist,abs(angdiffDeg([upperLimPerDistBin'; phaseAvgValuesPerDist'])),abs(angdiffDeg([lowerLimPerDistBin'; phaseAvgValuesPerDist'])),'k')
        %plot(currFieldDistInFieldFracPerSpike,currFieldSpikePhases,'k.')
        %}
        
        drawnow
        
         if(getBinAvgs)
              xlim([0 1])
            ylim([0 360])
         else
              xlim([0.1 1])
             ylim([0.1 1])
             axis square
             set(gca,'xscale','log')
                set(gca,'yscale','log')
                %close all
                
         end

         if(fi==100)
             disp('')
         end
     
    end
    
    
end
%%
save('allDistTimePhaseTriples1227FieldsWithMaxAcc.mat','minSpeed','maxSpeed','maxAccMPerSecSquared','usedFieldCount','minAccMPerSecSquared','allDistTimePhaseTriples','totalFieldCount','-v7.3')

%%

close all
input.data3cols=allDistTimePhaseTriples;

minNumPts=50;

numBins=50;
useLogBins=1;
useLogBins=0;
%useLogBins=0;

useLogTimeBins=0;

input.xRad=0.01;
input.yRad=0.01;

maxNegExp=1;
   
if(useLogBins)


    input.xEdges=10.^(-linspace(0,maxNegExp,numBins+1));
    input.yEdges=10.^(-linspace(0,maxNegExp,numBins+1));
else

    input.xEdges=linspace(0.05,1,numBins+1);
    input.yEdges=linspace(0.05,1,numBins+1);
end

if(useLogTimeBins)
   
      input.yEdges=10.^(-linspace(0,maxNegExp,numBins+1));
        input.xEdges=linspace(1,10^-maxNegExp,numBins+1);
end
input.desiredStat='circMean';

[xCenters,yCenters,allDistTimePhaseMeansSmooth,dataCountPerBin] =makeHeatMapOf3rdVarUniv(input);

figure;

allDistTimePhaseMeansSmoothDisp=allDistTimePhaseMeansSmooth;

allDistTimePhaseMeansSmoothDisp(allDistTimePhaseMeansSmooth<minNumPts)=NaN;
omarPcolor(xCenters,yCenters,allDistTimePhaseMeansSmoothDisp)
colormap(jet)
cb=colorbar

xlabel('Distance in field (frac)')

ylabel('Time in field (frac)')

ylabel(cb,'Circ mean phase (deg)')

title(sprintf('All dist-time-phase triplets, HC3 (n=%d fields)',totalFieldCount))

axis square
maxFig
setFigFontTo(16)


if(useLogBins)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end
if(useLogTimeBins)
    set(gca,'yscale','log')
end

if(useLogTimeBins)
    saveas(gcf,'hc3SpaceLogTimePhase.png')
elseif(useLogBins)
    saveas(gcf,'hc3LogSpaceLogTimePhase.png')
else
    saveas(gcf,'hc3SpaceTimePhase.png')
end





%omarPcolor(xCenters,yCenters,dataCountPerBin)
%colormap(parula)
%colorbar

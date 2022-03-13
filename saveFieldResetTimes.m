clear all; close all; clc

%codeVersionID='4-Jan-2021_6-12Hz_OnlyPrenatalAlcExposure';numFiles=61;
%codeVersionID='6-12Hz_OnlyPrenatalAlcExposure';

%codeVersionID='4-Jan-2021_6-12Hz_NoPrenatalAlcExposure';numFiles=23;
codeVersionID='6-12Hz_NoPAE';

if(~exist(sprintf('accumulatedPhaseData_%s.mat',codeVersionID),'file'))
    
imageDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong/tifs';
dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong';
%dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong/actuallyPAE';


saveDir='/Users/tibinjohn/thetaSeq/code/PhaseAndDistVsPhaseAndSpeedChange';

dataFilePaths=getFilePathsRegex(dataDir,'*mat');
imageFilePaths=getFilePathsRegex(imageDir,'*tif');

touchDir(saveDir)

masterDistPerTime=[];
masterLapTimePerTime=[];
masterSpikeRatePerTime=[];
masterSpikePhasePerTime=[];
masterSpeedPerTime=[];

masterInFieldDelXTP=[];
masterInFieldDelXTC=[];
masterInFieldDelXTS=[];



normDistWithinField=1;
normTimeWithinField=1;
normalizeByTrueLapAvgPeakRate=1;

%normDistWithinField=0;
%normTimeWithinField=0;
distRelativeToFieldStart=0;

normEachCell=1;

if(normTimeWithinField)
    %maxDispTime=1.2;
     maxDispTime=1;
    %maxDispTime=2;
else
    maxDispTime=3;

end
%maxDispTime=5;
numFiles=length(dataFilePaths);

for fi=1:numFiles
    currFilePath=dataFilePaths{fi};
    
    if(~contains(currFilePath,'LEM3206_S2019072221142_cellNum10_dir1'))
        
    %if(~contains(currFilePath,'LEM3246_S2019070914314_cellNum10_dir1'))
    %if(~contains(currFilePath,'LEM3246_S2019062916424_cellNum10_dir1'))
    %if(~contains(currFilePath,'LEM3246_S2019072917264_cellNum15_dir2'))
    %if(~contains(currFilePath,'LEM3206_S2019072221142_cellNum36_dir2'))
    
        %continue
    end
    
   
    currImgPath=imageFilePaths{fi};
    currData=load(currFilePath)
    
    currStartCm=currData.manualFieldStartCm;
    currEndCm=currData.manualFieldEndCm;
    %currStartCm=currData.manualFieldStartCmAdjusted;
    %currEndCm=currData.manualFieldEndCmAdjusted;
    manualPrecessionStartPhaseDeg=currData.manualPrecessionStartPhaseDeg;
    
    %cellIDstr=currData.spaceTimePhaseInfo.cellIDstr;
    fileName=getFileNameFromPath(currFilePath);
    cellIDstr=fileName(1:(end-28));
    
    if(currEndCm<=115)
        %continue
    end
   [placeCellDistPerTime,placeCellTimeInLapPerTime,placeCellRatePerTime,comDist,comTime]=addSpeedChangeToDataStruct(currFilePath);
   %continue
   newData=load(currFilePath);
   
   flipManualPlaceBounds=newData.youShouldFlipManualPlaceBounds; %always 0 now
   trackLengthCm=newData.trackLengthCm;
   
   if(flipManualPlaceBounds)
       currStartCm=trackLengthCm-currStartCm;
       currEndCm=trackLengthCm-currEndCm;
   end
   
   placeCellSpeedPerTime=newData.placeCellSpeedPerTime;
   
   avgPeakRateAcrossLaps=newData.avgPeakRateAcrossLaps;
   
   placeCellPhasePerTime=newData.placeCellPhasePerTime;
   placeCellThetaAmpPerTime=newData.placeCellThetaAmpPerTime;
   placeCellThetaDurPerTime=newData.placeCellThetaDurPerTime;
   
   placeCellDistPerTime=newData.placeCellFieldFracPerTime;
   
   placeCellTimeInLapPerTime=newData.placeCellTimeInFieldPerTime;
   
   avgWaveform=newData.avgWaveform;
   
   inFieldDelXTP=newData.inFieldDelXTP;
   
   inFieldDelXTP(:,3)=inFieldDelXTP(:,3)+180;
   inputXTP.data3cols=inFieldDelXTP;
    
   %inputXTP.desiredStat='mean';
   inputXTP.desiredStat='circMean';
    %inputXTP.xRad=1; %cm
    %inputXTP.xRad=1; %cm
     %inputXTP.xRad=1.5; %cm
     %inputXTP.xRad=1.5; %cm
    %inputXTP.yRad=0.015; %sec
   
    %inputXTP.xRad=1; %cm
    %inputXTP.yRad=0.01; %sec
    
    
    if(normEachCell)
        %minDistXTP=0.05;%fraction
        %maxDistXTP=0.65;%fraction
        minDistXTP=0;%fraction
        maxDistXTP=0.6;%fraction
        inputXTP.xRad=(maxDistXTP-minDistXTP)/15; %fraction
        

         %minTimeXTP=0.05;%sec
        %maxTimeXTP=0.65; %sec
          minTimeXTP=0;%sec
        maxTimeXTP=0.6; %sec
        
        
        %minTimeXTP=0.15;%sec
        %maxTimeXTP=0.35; %sec
        inputXTP.yRad=(maxTimeXTP-minTimeXTP)/15;%fraction
    else
        minDistXTP=0;%cm
        maxDistXTP=25;%cm
        inputXTP.xRad=1*2; %cm
        
        minTimeXTP=0.15;%sec
        maxTimeXTP=0.35; %sec
        inputXTP.yRad=0.01*2; %sec
    end
    
  
    
    inputXTP.xEdges=linspace(minDistXTP,maxDistXTP,75+1);
     inputXTP.yEdges=linspace(minTimeXTP,maxTimeXTP,75+1);
     
   [delDistBins,delTimeBins,delPhaseSmooth,dataCountPerBinXTP]=makeHeatMapOf3rdVarUniv(inputXTP);
   
   inFieldDelXTC=newData.inFieldDelXTC;
   inputXTP.data3cols=inFieldDelXTC;
      [delDistBinsCycleDurs,delTimeBinsCycleDurs,delCycleDursSmooth,dataCountPerBinXTC]=makeHeatMapOf3rdVarUniv(inputXTP);

   inFieldDelXTS=newData.inFieldDelXTS;
    inputXTP.data3cols=inFieldDelXTS;
       [delDistBinsSpeed,delTimeBinsSpeed,delSpeedSmooth,dataCountPerBinXTS]=makeHeatMapOf3rdVarUniv(inputXTP);

   
   delPhaseSmoothed=nanGaussSmooth(delPhaseSmooth);
   delCycleDursSmoothed=nanGaussSmooth(delCycleDursSmooth);
   delSpeedSmoothed=nanGaussSmooth(delSpeedSmooth);
   dataCountPerBinXTPsmoothed=nanGaussSmooth(dataCountPerBinXTP);
   %delPhaseSmoothed=delPhaseSmooth;
   fH=figure
  
   %subplot(2,2,1)
   uberTitle('test')
         ax=gca;
         
   plot(newData.allElapsedDistsInFieldCm/newData.meanFieldWidthCm...
           ,newData.allElapsedTimesInFieldSec/newData.meanFieldTimeWidthSec,'k-','LineWidth',3);
       hold on
       
       plotSpacetimePredictedPhaseCone
       
       xlabel('Distance in field (fields)')
       %ylabel('Speed in field (cm/s)')
       ylabel('Time in field (fields)')
         xlim([0 1])
       %ylim([5 max(newData.dispSpeedsInField)])
       ylim([0 1])
       
       ax.Position=[0 0.6 0.35 0.35];
       %ax.Position=[-0.1 0.25 0.55 0.55];
       
       %title('Distance-Speed profile in field')
       title('Spacetime path in field and phase cone model')
       
   uberTitle('test')
         ax=gca;
   omarPcolor(delDistBins,delTimeBins,(-(delPhaseSmoothed-180)')/360,fH)
        %xlabel('Distance along track (cm)')
        if(normEachCell)
             xlabel('\Deltax per cycle (fraction of field)')
             ylabel('\Deltat per cycle (fraction of field)')
        else
            xlabel('\Deltax (cm)')
             ylabel('\Deltat (sec)')
        end
        colormap(jet)
        cb=colorbar
        %ylabel(cb,'Mean \Deltap precession (degrees)')
        ylabel(cb,'Mean \Deltap precession per cycle (fraction of theta cycle)')
       %ylim([0 maxDispTime])
      %caxis([0 360])
      %caxis([0 180])
       %caxis([0 100])
        caxis([0 0.6])
        %caxis([0 1/7])
      xlim([minDistXTP maxDistXTP])
      %ylim([0 0.25])
      ylim([minTimeXTP maxTimeXTP])
       daspect([1 1 1])
       hold on
   
          title(sprintf('Average theta phase change in delta-cycle spacetime'))
           ax.Position=[0 0.1 0.35 0.35];
   %}
   
          %{
       subplot(2,2,2)
          omarPcolor(delDistBins,delTimeBins,dataCountPerBinXTPsmoothed',fH)
        %xlabel('Distance along track (cm)')
        if(normEachCell)
             xlabel('\Deltax (fraction of field)')
             ylabel('\Deltat (fraction of field)')
        else
            xlabel('\Deltax (cm)')
             ylabel('\Deltat (sec)')
        end
        colormap(gca,parula)
        cb=colorbar
        ylabel(cb,'Count')

      xlim([minDistXTP maxDistXTP])
      ylim([minTimeXTP maxTimeXTP])
       daspect([1 1 1])
       hold on
       %}
       
       %{
       subplot(2,2,3)
          omarPcolor(delDistBinsCycleDurs,delTimeBinsCycleDurs,1./delCycleDursSmoothed',fH)
        %xlabel('Distance along track (cm)')
        if(normEachCell)
             xlabel('\Deltax (fraction of field)')
             ylabel('\Deltat (fraction of field)')
        else
            xlabel('\Deltax (cm)')
             ylabel('\Deltat (sec)')
        end
        colormap(gca,parula)
        cb=colorbar
        caxis([5 10])
        ylabel(cb,'LFP theta frequency (Hz)')

      xlim([minDistXTP maxDistXTP])
      ylim([minTimeXTP maxTimeXTP])
       daspect([1 1 1])
       %}
       
       %{
       subplot(2,2,4)
          omarPcolor(delDistBinsSpeed,delTimeBinsSpeed,delSpeedSmoothed',fH)
        %xlabel('Distance along track (cm)')
        if(normEachCell)
             xlabel('\Deltax (fraction of field)')
             ylabel('\Deltat (fraction of field)')
        else
            xlabel('\Deltax (cm)')
             ylabel('\Deltat (sec)')
        end
        colormap(gca,parula)
        cb=colorbar
        caxis([0 50])
        ylabel(cb,'Speed (cm/s)')

      xlim([minDistXTP maxDistXTP])
      ylim([minTimeXTP maxTimeXTP])
       daspect([1 1 1])
       %}
       %plot([0 1],[1 0],'k--','LineWidth',3)
       
       
        uberTitle('test')
         ax=gca;
  
       plot(newData.allElapsedDistsInFieldCm,newData.dispSpeedsInField,'k-');
       
       
       %xlabel('Distance in field (fields)')
       ylabel('Speed in field (cm/s)')
       %ylabel('Time in field (fields)')
         %xlim([0 1])
       %ylim([5 max(newData.dispSpeedsInField)])
       %ylim([0 1])
       %ax.Position=[0.5 0.55 0.10 0.15];
             ax.Position=[0.675 0.8 0.10 0.15];
       title('Speed vs distance in field')
       %title('Spacetime path in field')
       commonDistLim=[0 max(newData.allElapsedDistsInFieldCm)];
       
       commonSpeedLim=[min(newData.dispSpeedsInField) max(newData.dispSpeedsInField)];
       ylim(commonSpeedLim)
       xlim([0 max(newData.allElapsedDistsInFieldCm)])
       
       uberTitle('test')
         ax=gca;
  
       %plot(newData.allElapsedTimesInFieldSec,newData.dispSpeedsInField,'k-');
         plot(newData.allElapsedDistsInFieldCm,newData.speedChangesPerPosDuringField,'k.');
       %xlabel('Time in field (sec)')
       xlabel('Distance in field (cm)')
       %ylabel('Speed in field (cm/s)')
       %ylabel('Speed change in field (cm/s/cm)')
      ylabel('Acc in field (cm/s^2)')
      
       %xlim([-Inf newData.meanFieldTimeWidthSec])
       %ylim(commonSpeedLim)
       xlim([0 max(newData.allElapsedDistsInFieldCm)])
       ylim([-3 3])
       
         
          ax.Position=[0.675 0.55 0.10 0.15];
       %title('Time-Speed profile in field')
        title('Acc. vs distance in field')
        
        
        
            uberTitle('test')
         ax=gca;
  
       %plot(newData.allElapsedTimesInFieldSec,newData.dispSpeedsInField,'k-');
         plot(newData.allElapsedDistsInFieldCm,newData.jerkDuringField,'k.');
       %xlabel('Time in field (sec)')
       xlabel('Distance in field (cm)')
      ylabel('Jerk in field (cm/s^2)')

       xlim([0 max(newData.allElapsedDistsInFieldCm)])
       ylim([-1 1])
       
         
         ax.Position=[0.675 0.3 0.10 0.15];
         %ax.Position=[0.5 0.3 0.10 0.15];
       %title('Time-Speed profile in field')
        title('Jerk vs distance in field')
        
       
         uberTitle('test')
         ax=gca;
  
        cycleChangeIdxes=newData.cycleChangeIdxesInField==1;

       plot(newData.allElapsedDistsInFieldCm(cycleChangeIdxes),newData.dispPhasesInField(cycleChangeIdxes)/360,'k.');
       xlabel('Distance in field (cm)')
       %ylabel('Theta phase in field (deg)')
       ylabel('Theta phase in field (cycle frac)')
         xlim(commonDistLim)
        %ylim([0 360])
        ylim([0 1])
       ax.Position=[0.5 0.8 0.10 0.15];
       title('Distance-Phase profile in field')
       
           cycleStepSizePerTime=inFieldDelXTP(:,1)./inFieldDelXTP(:,2);
        xChangePerTimeInField=inFieldDelXTP(:,1);
        tChangePerTimeInField=inFieldDelXTP(:,2);
        phaseStepSizePerTime=inFieldDelXTP(:,3)-180;
        
         %cycleStepSizePerTime=newData.twoCycleSpaceDiffPerTime./newData.twoCycleTimeDiffPerTime;
         dispCycIdxes=~isnan(cycleStepSizePerTime);
         dispXIdxes=~isnan(xChangePerTimeInField);
         dispTIdxes=~isnan(tChangePerTimeInField);
         dispPidxes=~isnan(phaseStepSizePerTime);
       
       uberTitle('test')
         ax=gca;
         
         %{
         %plot(inFieldDelXTS(dispXIdxes,3),phaseStepSizePerTime(dispXIdxes)/360./xChangePerTimeInField(dispXIdxes),'k.');
         %plot( ((xChangePerTimeInField(dispXIdxes)).*(newData.speedChangesPerPosDuringField(dispXIdxes))),...
         %    phaseStepSizePerTime(dispXIdxes)/360,'k.');
           plot( xChangePerTimeInField(dispXIdxes),...
             phaseStepSizePerTime(dispXIdxes)/360,'k.');
         xlabel('\DeltaX per cycle')
        ylabel('\DeltaP per cycle')
        
      %xlim(commonSpeedLim)
          %xlim([-1.5 1.5])
          xlim([0 0.6])
          
          commonDeltaPLim=[-0.6 0];
      ylim(commonDeltaPLim)
      hold on
      plot(xlim,-ylim,'k--')
      title('\Deltaphase vs \DeltaX in field')
      %}
 
         
            plot(newData.dispSpeedsInField(dispXIdxes),phaseStepSizePerTime(dispXIdxes)/360./xChangePerTimeInField(dispXIdxes),'k.');
             ylabel('\DeltaP/\DeltaX per cycle')
             xlabel('Speed per cycle (cm/s)')
      xlim(commonSpeedLim)
      ylim([-1.5 0]) 
               title('\Deltaphase/\DeltaX vs speed in field')

      ax.Position=[0.85 0.8 0.10 0.15];
      %ax.Position=[0.5 0.55 0.10 0.15];
       
  
         %{
       plot(newData.allElapsedTimesInFieldSec(cycleChangeIdxes),newData.dispPhasesInField(cycleChangeIdxes),'k.');
       xlabel('Time in field (sec)')
       ylabel('Theta phase (deg)')
      
          xlim([-Inf newData.meanFieldTimeWidthSec])
       ylim([0 360])
            ax.Position=[0.85 0.8 0.10 0.15];
       title('Time-Phase profile in field')
         %}
       
       
       uberTitle('test')
         ax=gca;
  
         %cycleStepSizePerTime=inFieldDelXTC(:,1)./inFieldDelXTC(:,2);
         %%cycleStepSizePerTime=newData.twoCycleSpaceDiffPerTime./newData.twoCycleTimeDiffPerTime;
         %dispCycIdxes=~isnan(cycleStepSizePerTime);
        %plot(cycleStepSizePerTime(dispCycIdxes),1./newData.dispCycDursInField(dispCycIdxes),'k.');
        
    
         
        %plot(cycleStepSizePerTime(dispCycIdxes),phaseStepSizePerTime(dispCycIdxes)/360,'k.');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %deltaP deltaX ratio plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(newData.allElapsedDistsInFieldCm(dispXIdxes)...
            ,phaseStepSizePerTime(dispXIdxes)/360./xChangePerTimeInField(dispXIdxes),'k.');
     
         %plot(1./newData.dispCycDursInField(dispPidxes),phaseStepSizePerTime(dispPidxes)/360,'k.');

       %plot(cycleStepSizePerTime(cycleChangeIdxes),1./newData.dispCycDursInField(cycleChangeIdxes),'k.');
       %xlabel('\DeltaX/\DeltaT per cycle')
       %xlabel('\DeltaX per cycle (field frac)')
       %ylabel('\Theta freq (Hz)')
        %ylabel('\DeltaP(fraction of cycle)')
        %{
         plot( newData.dispSpeedsInField(dispXIdxes)/newData.meanFieldWidthCm,...
             phaseStepSizePerTime(dispXIdxes)/360,'k.');
         xlabel('Speed per cycle (fields/s)')
         %xlim(commonSpeedLim/newData.meanFieldWidthCm)
          xlim([0.1 1.1])
        ylabel('\DeltaP per cycle')
        ylim(commonDeltaPLim)
         title('\Deltaphase vs Speed in field')
        %}
        
         
         ylabel('\DeltaP/\DeltaX per cycle')
          xlabel('Distance in field (cm)')
        xlim([0 max(newData.allElapsedDistsInFieldCm)])
      
          %xlim([-Inf newData.meanFieldTimeWidthSec])
          %xlim([0 0.6])
    %xlim([prctile(xChangePerTimeInField,1) prctile(xChangePerTimeInField,99)])
        ylim([-1.5 0])
       %ylim([5 12])
       %ylim([-180 0])
        %ylim([-180 0]/360)
        %ylim([-60 60])
           %ylim([-60 0])
   
             ax.Position=[0.85 0.8 0.10 0.15];
             ax.Position=[0.5 0.55 0.10 0.15];
            %xlim([5 10])
       %title('Cycle step size vs \Theta Freq in field')
       title('\Deltaphase/\DeltaX vs dist in field')
      
       
          uberTitle('test')
         ax=gca;
  
       %plot(newData.allElapsedTimesInFieldSec(cycleChangeIdxes),newData.dispSpeedsInField(cycleChangeIdxes).*newData.dispCycDursInField(cycleChangeIdxes)/newData.meanFieldWidthCm,'k.');
      %plot(cycleStepSizePerTime(dispCycIdxes),newData.dispSpeedsInField(dispCycIdxes)/newData.meanFieldWidthCm,'k.');
     %plot(tChangePerTimeInField(dispTIdxes),phaseStepSizePerTime(dispTIdxes)/360,'k.');
   %{
     plot( newData.speedChangesPerPosDuringField(dispXIdxes)/newData.meanFieldWidthCm,...
             phaseStepSizePerTime(dispXIdxes)/360,'k.');
         xlabel('Acc. per cycle (fields/s^2)')
         %xlim([-2.5 2.5]/newData.meanFieldWidthCm)
         xlim([-0.05 0.05])
        ylabel('\DeltaP per cycle')
        ylim(commonDeltaPLim)
     %}
        
     plot(newData.speedChangesPerPosDuringField(dispXIdxes),phaseStepSizePerTime(dispXIdxes)/360./xChangePerTimeInField(dispXIdxes),'k.');
             ylabel('\DeltaP/\DeltaX per cycle')
             xlabel('Acc. per cycle (cm/s^2)')
      xlim([-3 3])
      ylim([-1.5 0]) 
               title('\Deltaphase/\DeltaX vs acc. in field')
         
       %ylabel('\Theta freq (Hz)')
        %ylabel('\DeltaP (cycle frac))')
   
       %xlabel('Time in field (sec)')
              %xlabel('\DeltaX/\DeltaT per cycle')
     
       %ylabel('Speed/\Thetafreq (fields/cycle)')
      %ylabel('Speed (fields/sec)')
       %ylabel('Speed (fields/sec)')
         %ylabel('\DeltaP(fraction of cycle)')
          %xlim([-Inf newData.meanFieldTimeWidthSec])
           %xlim([min(tChangePerTimeInField) max(tChangePerTimeInField)])
          %xlim([prctile(tChangePerTimeInField,1) prctile(tChangePerTimeInField,99)])
          %xlim([5 max(newData.dispSpeedsInField)])
             %xlim([0 0.6])
       %ylim([5 12])
       %ylim([0 1])
       %ylim([-180 0]/360)
       
      
         ax.Position=[0.85 0.55 0.10 0.15];
       %title('(Speed/\ThetaFreq) in field')
         %title('Cycle step size vs Speed in field')
         %title('\DeltaT vs \Deltaphase in field')

         %title('\Deltaphase vs Acc. in field')
       
           uberTitle('test')
         ax=gca;
  
         
            plot(newData.jerkDuringField(dispXIdxes),phaseStepSizePerTime(dispXIdxes)/360./xChangePerTimeInField(dispXIdxes),'k.');
             ylabel('\DeltaP/\DeltaX per cycle')
             xlabel('Jerk per cycle (cm/s^3)')
      xlim([-1 1])
      ylim([-1.5 0]) 
               title('\Deltaphase/\DeltaX vs jerk in field')
               ax.Position=[0.85 0.3 0.10 0.15];
               
               
         uberTitle('test')
         ax=gca;
  
         
         xVals=newData.allElapsedDistsInFieldCm/newData.meanFieldWidthCm;
         yVals=newData.allElapsedTimesInFieldSec/newData.meanFieldTimeWidthSec;
         
         xyPairs=[xVals(:) yVals(:)];
         
         [Zq] = plotSpacetimePredictedPhaseCone(xyPairs);
            plot(newData.allElapsedDistsInFieldCm/newData.meanFieldWidthCm,Zq,'k.');
             ylabel('Predicted phase (cycle frac)')
             xlabel('Distance in field (fields)')
      xlim([0 1])
      ylim([0 1]) 
               title('Spacetime cone-predicted phase in field')
               ax.Position=[0.5 0.3 0.10 0.15];
               
        
       
       maxFig
       setFigFontTo(18) 
   
       print('-r75',sprintf('deltaSpacetimeVsPhaseChange_%s_%s',cellIDstr,codeVersionID),'-dpng')
       close all
   
   
      minSpeedDispSpikes=2; %cm/s
   minThetaAmpDispSpikes=2; %LFP Z
   minSpeedIdxes=placeCellSpeedPerTime>minSpeedDispSpikes;
   minThetaIdxes=placeCellThetaAmpPerTime>minThetaAmpDispSpikes;
   
   minSpeedIdxes=minSpeedIdxes(:);
   minThetaIdxes=minThetaIdxes(:);
   
   goodIdxesSpikes=minSpeedIdxes & minThetaIdxes;

   maxLapTime=newData.maxLapTime;
   
   

   inputStruct.xRad=3; %cm
   %inputStruct.xRad=1; %cm
   %inputStruct.yRad=0.1 %sec
   
   %too big of squares...
    inputStruct.xRad=3; %cm
    inputStruct.yRad=0.3 %sec
    
     
    inputStruct.xRad=1; %cm
    inputStruct.yRad=0.1 %sec
    
    %inputStruct.xRad=1; %cm
   %inputStruct.yRad=0.15 %sec
    
       inputStruct.xRad=1; %cm
    inputStruct.yRad=0.3 %sec
    
    
    
       %inputStruct.xRad=0.5; %cm
   % inputStruct.yRad=0.3 %sec
    
    %too smooth
    %inputStruct.xRad=2; %cm
    %inputStruct.yRad=0.2 %sec

    inputStruct.xRad=1; %cm
    inputStruct.xEdges=linspace(0,trackLengthCm,75+1);
   if(normDistWithinField)
    inputStruct.xRad=0.01; %field frac
   inputStruct.xEdges=linspace(0,1*1,75/1+1);
 
   end
   %inputStruct.yRad=0.3 %sec
   inputStruct.yRad=0.1 %sec
    inputStruct.yEdges=linspace(0,maxLapTime,75+1);
    
    if(normTimeWithinField)
        %inputStruct.yRad=0.1/3; %field frac
        %inputStruct.yEdges=linspace(0,1,75+1);
    
    %inputStruct.yRad=0.01; %field frac
    inputStruct.yRad=0.1/3; %field frac
    %inputStruct.yEdges=linspace(0,1.2,(75/1.2+1));
    inputStruct.yEdges=linspace(0,1,(75+1));
        
    end
        
   
   if(distRelativeToFieldStart)
       inputStruct.xRad=1; %field frac
        inputStruct.yRad=0.3 %sec
       inputStruct.xEdges=linspace(0,trackLengthCm,75+1);
       inputStruct.yEdges=linspace(0,maxLapTime,75+1);
   end
   
   %inputStruct.xEdges=linspace(0,trackLengthCm,50+1);
   %inputStruct.yEdges=linspace(0,maxLapTime,50+1);
   %inputStruct.xEdges=linspace(0,trackLengthCm,100+1);
   %inputStruct.yEdges=linspace(0,maxLapTime,100+1);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %spike rate stats spacetime heat map
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if(normalizeByTrueLapAvgPeakRate)
        data3cols=[placeCellDistPerTime(:),placeCellTimeInLapPerTime(:),placeCellRatePerTime(:)/avgPeakRateAcrossLaps];
   else
       data3cols=[placeCellDistPerTime(:),placeCellTimeInLapPerTime(:),placeCellRatePerTime(:)];
   end
   inputStruct.data3cols=data3cols;
   inputStruct.desiredStat='mean';
   [distBins,timeBins,meanRateSmooth,dataCountPerBin]=makeHeatMapOf3rdVarUniv(inputStruct);
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %spike phase stats spacetime heat map
   %include only non-noise level theta-associated spikes
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %data3cols=[placeCellDistPerTime(:),placeCellTimeInLapPerTime(:),placeCellPhasePerTime(:)];
   
   highThetaDists=placeCellDistPerTime(minThetaIdxes);
   highThetaTimes=placeCellTimeInLapPerTime(minThetaIdxes);
   highThetaPhases=placeCellPhasePerTime(minThetaIdxes);
   data3cols=[highThetaDists(:), highThetaTimes(:),highThetaPhases(:)];
   
   
   inputStruct.data3cols=data3cols;
   inputStruct.desiredStat='circMean';
   [distBinsPhase,timeBinsPhase,meanPhaseSmooth,dataCountPerBinPhase]=makeHeatMapOf3rdVarUniv(inputStruct);
   
   inputStruct.desiredStat='circR';
   [distBinsPhase,timeBinsPhase,meanPhaseRSmooth,dataCountPerBinPhaseR]=makeHeatMapOf3rdVarUniv(inputStruct);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %theta stats spacetime heat map
   %include only non-noise level theta
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   highThetaAmps=placeCellThetaAmpPerTime(minThetaIdxes);

   data3cols=[highThetaDists(:), highThetaTimes(:),highThetaAmps(:)];

   inputStruct.data3cols=data3cols;
   inputStruct.desiredStat='mean';
   [distBinsTheta,timeBinsTheta,meanCycleAmp,dataCountPerBinTheta]=makeHeatMapOf3rdVarUniv(inputStruct);
   
   highThetaDur=placeCellThetaDurPerTime(minThetaIdxes);
   data3cols=[highThetaDists(:), highThetaTimes(:),highThetaDur(:)];
   
    inputStruct.data3cols=data3cols;
   inputStruct.desiredStat='mean';
   [distBinsTheta,timeBinsTheta,meanCycleDur,dataCountPerBinTheta]=makeHeatMapOf3rdVarUniv(inputStruct);
   
   %{
   scatterDists=placeCellDistPerTime(goodIdxesSpikes);
   scatterTimes=placeCellTimeInLapPerTime(goodIdxesSpikes);
   scatterPhases=placeCellPhasePerTime(goodIdxesSpikes);
   
    figure
   scatter(scatterDists,scatterTimes,25,scatterPhases,'filled')
   colormap(jet)
   colorbar
   set(gca,'Color','k')
   disp('')
   %}
   
   %%
   
   

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %speed stats spacetime heat map
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   data3cols=[placeCellDistPerTime(:),placeCellTimeInLapPerTime(:),placeCellSpeedPerTime(:)];
   inputStruct.data3cols=data3cols;
   
   inputStruct.desiredStat='mean';
   [distBinsSpeed,timeBinsSpeed,meanSpeedSmooth,dataCountPerBinSpeed]=makeHeatMapOf3rdVarUniv(inputStruct);
  
   inputStruct.desiredStat='SEM';
   %[distBinsSpeed,timeBinsSpeed,semSpeedSmooth,dataCountPerBinSpeed]=makeHeatMapOf3rdVarUniv(inputStruct);

   %{
     masterNumBins=50;
   masterNumBins=33;
   
   [distBins,timeBins,meanRateSmooth,dataCountPerBin] = makeHeatMapAvgOf3rdVar(placeCellDistPerTime,placeCellTimeInLapPerTime,placeCellRatePerTime,masterNumBins);
   [distBinsPhase,timeBinsPhase,meanPhaseSmooth,dataCountPerBinPhase] = makeHeatMapCircAvgOf3rdVar(placeCellDistPerTime,placeCellTimeInLapPerTime,placeCellPhasePerTime,masterNumBins);
   [distBinsPhase,timeBinsPhase,meanPhaseRSmooth,dataCountPerBinPhaseR] = makeHeatMapCircROf3rdVar(placeCellDistPerTime,placeCellTimeInLapPerTime,placeCellPhasePerTime,masterNumBins);
   [distBinsSpeed,timeBinsSpeed,meanSpeedSmooth,dataCountPerBinSpeed] = makeHeatMapAvgOf3rdVar(placeCellDistPerTime,placeCellTimeInLapPerTime,placeCellSpeedPerTime,masterNumBins);
    %}
   
   fH=figure
   
    
   %{
     subplot(2,2,2)
        omarPcolor(distBins+comDist,timeBins+comTime,dataCountPerBin',fH)
        xlabel('Distance along track (cm)')
        ylabel('Time in true lap (sec)')
        colormap(parula)
        cb=colorbar
        ylabel(cb,'data point count')
           daspect([1 1 1])
           
      %}
         
           
    %subplot(2,2,[1 2])
     subplot(3,2,[1])
    omarPcolor(distBins+comDist,timeBins+comTime,meanRateSmooth',fH)

         %xlabel('Distance along track (cm)')
         xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
    colormap(gca,parula)
    cb=colorbar
    if(normalizeByTrueLapAvgPeakRate)
        ylabel(cb,'Avg firing rate (fraction of peak)')
    else
    ylabel(cb,'Avg firing rate (Hz)')
    end
    %title(sprintf('Direction %d',di))
    
    title('Place cell firing rate')
    ylim([0 maxDispTime])
 

      commonXlim=xlim;
      commonYlim=ylim;
      
        
    hold on
     plot([currStartCm currStartCm],commonYlim,'k--','LineWidth',4)
        plot([currEndCm currEndCm],commonYlim,'k--','LineWidth',4)
        
      
           
    xlim(commonXlim)
    ylim(commonYlim)
     %set(gca, 'XTick', 0:5:120)
  %xlim([-30 30])
   % xlim([-50 50])
  %ylim([-1.5 1.5])
   % ylim([-3 3])
      if(normalizeByTrueLapAvgPeakRate)
           caxis([0 1])
          else
             caxis([0 30])
          end

    %ylim([0 maxDispTime])
   %daspect([1 1 1])
   
   
   subplot(3,2,[4])
     omarPcolor(distBinsTheta+comDist,timeBinsTheta+comTime,meanCycleAmp',fH)
        %xlabel('Distance along track (cm)')
        xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,parula)
        cb=colorbar
        ylabel(cb,'Mean theta amp (LFP Z)')
        caxis([2 5])
         xlim(commonXlim)
    ylim(commonYlim)
           %daspect([1 1 1])
    %}
    title('Theta amplitude')
   
        hold on
        
        plot([currStartCm currStartCm],commonYlim,'k--','LineWidth',4)
        plot([currEndCm currEndCm],commonYlim,'k--','LineWidth',4)
           
        
         %subplot(2,2,3)
         %subplot(2,2,[3 4])
          %subplot(4,2,[5 6])
          subplot(3,2,3)
         if(size(meanPhaseSmooth',1)<length(timeBinsPhase))
             timeBinsPhase=timeBinsPhase(1:size(meanPhaseSmooth',1));
         end
         
         if(size(meanPhaseSmooth',2)<length(distBinsPhase))
             distBinsPhase=distBinsPhase(1:size(meanPhaseSmooth',2));
         end
        omarPcolor(distBinsPhase+comDist,timeBinsPhase+comTime,meanPhaseSmooth',fH)
        
        hold on
        
        plot([currStartCm currStartCm],commonYlim,'k--','LineWidth',4)
        plot([currEndCm currEndCm],commonYlim,'k--','LineWidth',4)
     
   
     %xlabel('Distance along track (cm)')
     xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,jet)
        cb=colorbar
        ylabel(cb,'circ mean spike phase in theta')
        caxis([0 360])
         xlim(commonXlim)
    ylim(commonYlim)
    %ylim([0 maxDispTime])
    
    %set(gca, 'XTick', 0:5:120)
    
    title({'Place cell mean phase in local theta cycle'})
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %show consistency of phase in each space time bin
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %subplot(4,2,[7 8])
     subplot(3,2,[5])
     %omarPcolor(distBinsPhase+comDist,timeBinsPhase+comTime,meanPhaseKappaSmooth',fH)
      %omarPcolor(distBinsPhase+comDist,timeBinsPhase+comTime,meanPhaseCircSEMSmooth',fH)
    
      omarPcolor(distBinsPhase+comDist,timeBinsPhase+comTime,meanPhaseRSmooth',fH)
      
        
        hold on
        
        plot([currStartCm currStartCm],commonYlim,'k--','LineWidth',4)
        plot([currEndCm currEndCm],commonYlim,'k--','LineWidth',4)
     
   
     %xlabel('Distance along track (cm)')
     xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,parula)
        cb=colorbar
        ylabel(cb,'Resultant vector length spike phases')
        %caxis([0 360])
         xlim(commonXlim)
    ylim(commonYlim)
    %ylim([0 maxDispTime])
    caxis([0 1])
    %set(gca, 'XTick', 0:5:120)
    
    title('Place cell theta phase consistency')
     %daspect([1 1 1])
  %close all
  %xlim([-30 30])
   %  ylim([-2 2])
    %   xlim([-50 50])
     %ylim([-3 3])
  %caxis([0 30])

    
    %subplot(4,2,[3 4])
    subplot(3,2,2)
       omarPcolor(distBinsSpeed+comDist,timeBinsSpeed+comTime,meanSpeedSmooth',fH)
        %xlabel('Distance along track (cm)')
        xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,parula)
        cb=colorbar
        ylabel(cb,'Avg speed (cm/s)')
        caxis([0 40])
         xlim(commonXlim)
    ylim(commonYlim)
           %daspect([1 1 1])
    %}
  
   
        hold on
        
        plot([currStartCm currStartCm],commonYlim,'k--','LineWidth',4)
        plot([currEndCm currEndCm],commonYlim,'k--','LineWidth',4)
        title('Running Speed')
        
           subplot(3,2,6)
       omarPcolor(distBinsTheta+comDist,timeBinsTheta+comTime,1./meanCycleDur',fH)
        %xlabel('Distance along track (cm)')
        xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,parula)
        cb=colorbar
        ylabel(cb,'Avg theta freq (Hz)')
        %caxis([5 10])
        caxis([6 9])
         xlim(commonXlim)
    ylim(commonYlim)
           %daspect([1 1 1])
    %}
  
   
        hold on
        
        plot([currStartCm currStartCm],commonYlim,'k--','LineWidth',4)
        plot([currEndCm currEndCm],commonYlim,'k--','LineWidth',4)
        title('Theta frequency')
        

  
   setFigFontTo(14)
   
   maxFig
   uberTitle({'Spacetime behavior and hippocampal activity',removeUnderscores(cellIDstr),sprintf('Current field start: %d cm, field end: %d cm',currStartCm,currEndCm)})

   ax=gca;
   waveformTimeAxis=(1/32000)*(1:length(avgWaveform));
   pW=plot(waveformTimeAxis*1000,avgWaveform');
   xlabel('Time (msec)')
   xlim([-Inf Inf])
   ylim([-Inf Inf])
   ax.Position=[0.930 0.900 0.0600 0.0600];
   title('avg tetrode waveform')
   %saveas(gcf,sprintf('spacetimeRateField_%s.tif',cellIDstr))
   %print('-r100',sprintf('spacetimeRateField_%s',cellIDstr),'-dpng')
      uberTitle({'Spacetime behavior and hippocampal activity',removeUnderscores(cellIDstr),sprintf('Current field start: %d cm, field end: %d cm',currStartCm,currEndCm)})

    print('-r75',sprintf('spacetimeRateField_%s_%s',cellIDstr,codeVersionID),'-dpng')
   close all
   
    masterDistPerTime=[masterDistPerTime; placeCellDistPerTime(:)];
    masterLapTimePerTime=[masterLapTimePerTime; placeCellTimeInLapPerTime(:)];
    
     if(normalizeByTrueLapAvgPeakRate)
         masterSpikeRatePerTime=[masterSpikeRatePerTime;(placeCellRatePerTime(:)/avgPeakRateAcrossLaps)];
     else
        masterSpikeRatePerTime=[masterSpikeRatePerTime;placeCellRatePerTime(:)];
     end
    
    masterSpikePhasePerTime=[masterSpikePhasePerTime;placeCellPhasePerTime(:)];
    masterSpeedPerTime=[masterSpeedPerTime;placeCellSpeedPerTime(:)];
    
    masterInFieldDelXTP=[masterInFieldDelXTP;inFieldDelXTP];
    masterInFieldDelXTC=[masterInFieldDelXTC;inFieldDelXTC];
    masterInFieldDelXTS=[masterInFieldDelXTS;inFieldDelXTS];
   %addFieldResetTimesToDataStruct(currFilePath);


    fileSaveRootName=getFileNameFromPath(currFilePath);
    saveRootName=fileSaveRootName(1:(end-4));

    %[distBins,timeBins,meanRateSmooth] = makeHeatMapAvgOf3rdVar(masterDistPerTime,masterLapTimePerTime,masterSpikeRatePerTime,30);



end
%%
%[distBins,timeBins,meanRateSmooth,dataCountPerBin] = makeHeatMapAvgOf3rdVar(masterDistPerTime,masterLapTimePerTime,masterSpikeRatePerTime,300);



if(normDistWithinField)
    inputStruct.xRad=0.01; %field frac
    inputStruct.yRad=0.3 %sec
   %inputStruct.xEdges=linspace(0,1,75+1);
   %inputStruct.yEdges=linspace(0,maxLapTime,75+1);
    inputStruct.xRad=0.01; %field frac
   inputStruct.xEdges=linspace(0,1*1,75/1+1);
 
   
end

 if(normTimeWithinField)
        %inputStruct.yRad=0.1/3 %field frac
    %inputStruct.yEdges=linspace(0,1,75+1);
    %inputStruct.yRad=0.01; %field frac
    inputStruct.yRad=0.1/3; %field frac
    %inputStruct.yEdges=linspace(0,1.2,(75*1.2+1));
        
     inputStruct.yEdges=linspace(0,1,(75+1));
    end
   
if(distRelativeToFieldStart)
   inputStruct.xRad=1; %field frac
    inputStruct.yRad=0.3 %sec
   inputStruct.xEdges=linspace(0,trackLengthCm,75+1);
   inputStruct.yEdges=linspace(0,maxLapTime,75+1);
end
   

%inputStruct.xRad=0.5;


data3cols=[masterDistPerTime(:),masterLapTimePerTime(:),masterSpeedPerTime(:)];
   inputStruct.data3cols=data3cols;
   
   inputStruct.desiredStat='mean';
   [distBinsSpeed,timeBinsSpeed,meanSpeedSmooth,dataCountPerBinSpeed]=makeHeatMapOf3rdVarUniv(inputStruct);


  

   inputStruct.desiredStat='mean';
   inputStruct.data3cols=[masterDistPerTime(:), masterLapTimePerTime(:),masterSpikeRatePerTime(:)];
   
[distBins,timeBins,meanRateSmooth,dataCountPerBin] = makeHeatMapOf3rdVarUniv(inputStruct);

 inputStruct.desiredStat='circMean';
 inputStruct.data3cols=[masterDistPerTime(:), masterLapTimePerTime(:),masterSpikePhasePerTime(:)];
[distBinsPhase,timeBinsPhase,meanPhaseSmooth,dataCountPerBinPhase] = makeHeatMapOf3rdVarUniv(inputStruct);

 inputStruct.desiredStat='circR';
[distBinsPhase,timeBinsPhase,meanPhaseRSmooth,dataCountPerBinPhase] = makeHeatMapOf3rdVarUniv(inputStruct);

%{
goodCountIDs=dataCountPerBin>10;
goodCountIDs=dataCountPerBin>=3;
goodCountIDs=dataCountPerBin>=0;

meanRateSmooth(~goodCountIDs)=NaN;
%}

    save(sprintf('accumulatedPhaseData_%s.mat',codeVersionID),'-v7.3')
else
    
    load(sprintf('accumulatedPhaseData_%s.mat',codeVersionID))
%%
originalRate=meanRateSmooth;
originalPhase=meanPhaseSmooth;
originalPhaseR=meanPhaseRSmooth;
originalSpeed=meanSpeedSmooth;

meanRateSmoothed=nanGaussSmooth(originalRate);
meanPhaseSmoothed=nanGaussSmooth(originalPhase);
meanPhaseRSmoothed=nanGaussSmooth(originalPhaseR);
meanSpeedSmoothed=nanGaussSmooth(originalSpeed);


minSampleCount=30;
%minSampleCount=25;
%minSampleCount=20;

goodCountIDs=dataCountPerBin>=minSampleCount;
%goodCountIDs=dataCountPerBin>=3;
%goodCountIDs=dataCountPerBin>=0;
goodCountIDsPhase=dataCountPerBinPhase>=minSampleCount;
goodCountIDsSpeed=dataCountPerBinSpeed>=minSampleCount;

meanRateSmoothed(~goodCountIDs)=NaN;
meanPhaseSmoothed(~goodCountIDsPhase)=NaN;
meanPhaseRSmooth(~goodCountIDsPhase)=NaN;
meanSpeedSmoothed(~goodCountIDsSpeed)=NaN;


%meanRateSmoothed=originalRate;
%meanPhaseSmoothed=originalPhase;


   maxDispDist=1
   %maxDispTime=1.2
   maxDispTime=1
   %%
   figure
   subplot(1,2,1)
   [X,T]=meshgrid(distBinsPhase,timeBinsPhase);
   [C,hD]=contour(X,T,meanPhaseSmoothed'/360,35)
   hD.LineWidth=2;
   colormap(gca,jet)
   cb=colorbar
   caxis([0 1])
   xlabel('Distance in field (frac)')
   ylabel('Time in field (frac)')
   ylabel(cb,'Avg phase in spacetime (cycle frac)')
   title(sprintf('%d place cell grouped spacetime response',numFiles))
   
   hold on
plot([0 0.7],[0 1],'k--','LineWidth',3)
plot([0 1],[0 0.7],'k--','LineWidth',3)
%plot([0 1],[0 1],'k--','LineWidth',2)
box off
   daspect([1 1 1])
   subplot(1,2,2)
   plotSpacetimePredictedPhaseCone
   hold on 
   plot([0 0.7],[0 1],'k--','LineWidth',3)
plot([0 1],[0 0.7],'k--','LineWidth',3)

        uberTitle('test')
        ax=gca;
           % plot(0,0)
            plotSpacetimePredictedPhaseCone(0,'surface')
             
              % xlabel('Distance in field (frac)')
   %ylabel('Time in field (frac)')
   %zlabel('Theta phase (frac)')
   xticklabels({})
   yticklabels({})
   zticklabels({})
   
               ax.Position=[0.8 0.825 0.175 0.175];
               %axes(ax,'Color','none')
               view(-5.5,10.8)
               
uberTitle({'Theta phase as coding for height in spacetime cone','(speed-invariant hyperbolas)'})
maxFig
setFigFontTo(18)
   
saveas(gcf,sprintf('spacetimeConeDataAndModel%s.png',codeVersionID))
   %%
   
   
fH=figure
        subplot(2,2,1)
        omarPcolor(distBins,timeBins,meanRateSmoothed',fH)
        %xlabel('Distance along track (cm)')
        xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,jet)
        cb=colorbar
          if(normalizeByTrueLapAvgPeakRate)
            ylabel(cb,'Avg firing rate (fraction of peak)')
          else
            ylabel(cb,'Avg firing rate (Hz)')
          end
       ylim([0 maxDispTime])
       xlim([0 maxDispDist])
       
     if(normalizeByTrueLapAvgPeakRate)
           caxis([0 1])
           caxis([0 0.6])
          else
             caxis([0 30])
          end
       daspect([1 1 1])
   hold on
       
       plot([0 1],[1 0],'k--','LineWidth',3)
       
        title(sprintf('Average firing rate in spacetime (N>%d)',minSampleCount))
  
  %figure
   xlim([0 maxDispDist])
         subplot(2,2,2)
        %omarPcolor(distBins,timeBins,dataCountPerBin',fH)
        omarPcolor(distBinsSpeed,timeBinsSpeed,meanSpeedSmoothed',fH)
           %xlabel('Distance along track (cm)')
           xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        %xlabel('Distance from field COM (cm)')
        %ylabel('Time in true lap from field COM (sec)')
        colormap(gca,parula)
        cb=colorbar
        %ylabel(cb,'rate sample count')
         ylabel(cb,'Avg running speed (cm/s)')
   
           daspect([1 1 1])
        title(sprintf('Average running speedin spacetime (N>%d)',minSampleCount))
     ylim([0 maxDispTime])
  caxis([0 35])
  daspect([1 1 1])
     hold on
       
       plot([0 1],[1 0],'k--','LineWidth',3)
       
     xlim([0 maxDispDist])
    subplot(2,2,3)
        omarPcolor(distBinsPhase,timeBinsPhase,meanPhaseSmoothed',fH)
        %xlabel('Distance along track (cm)')
        xlabel('Distance in field (fraction)')
         if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(jet)
        cb=colorbar
        ylabel(cb,'Circ mean theta phase')
       ylim([0 maxDispTime])
      caxis([30 320])
      %caxis([100 300])
       daspect([1 1 1])
       hold on
       
       plot([0 1],[1 0],'k--','LineWidth',3)
       
       title(sprintf('Circ average theta phase in spacetime (N>%d)',minSampleCount))
       
           xlim([0 maxDispDist]) 
    subplot(2,2,4)
        omarPcolor(distBinsPhase,timeBinsPhase,meanPhaseRSmoothed',fH)
        %xlabel('Distance along track (cm)')
        xlabel('Distance in field (fraction)')
 if(normTimeWithinField)
             ylabel('Time in field (fraction)')
         else
             ylabel('Time (sec)')
         end
        colormap(gca,parula)
        cb=colorbar
        ylabel(cb,'phase resultant vector length')
       ylim([0 maxDispTime])
      caxis([0 0.6])
       daspect([1 1 1])
        title(sprintf('Theta phase consistency in spacetime (N>%d)',minSampleCount))
     hold on
       
     plot([0 1],[1 0],'k--','LineWidth',3)
       
       xlim([0 maxDispDist])
  
  
  
  
   setFigFontTo(18)
   
   uberTitle('place cell population summary')
   maxFig
   
       print('-r75',sprintf('spacetimeRateFieldSummary_%s',codeVersionID),'-dpng')
 %%      
       
       %%%%%%%%%%%%%%%%%
       %delta XTP
       %%%%%%%%%%%%%%%%%
       
       %masterInFieldDelXTP(:,3)=masterInFieldDelXTP(:,3)+180;
       inputXTP.data3cols=masterInFieldDelXTP;

   [delDistBins,delTimeBins,delPhaseSmooth,dataCountPerBinXTP]=makeHeatMapOf3rdVarUniv(inputXTP);
   
      
   inputXTP.data3cols=masterInFieldDelXTC;
      [delDistBinsCycleDurs,delTimeBinsCycleDurs,delCycleDursSmooth,dataCountPerBinXTC]=makeHeatMapOf3rdVarUniv(inputXTP);

   
    inputXTP.data3cols=masterInFieldDelXTS;
       [delDistBinsSpeed,delTimeBinsSpeed,delSpeedSmooth,dataCountPerBinXTS]=makeHeatMapOf3rdVarUniv(inputXTP);

   %%
   delPhaseSmoothed=nanGaussSmooth(delPhaseSmooth);
   delCycleDursSmoothed=nanGaussSmooth(delCycleDursSmooth);
   delSpeedSmoothed=nanGaussSmooth(delSpeedSmooth);
   
   %delPhaseSmoothed=delPhaseSmooth;
   
   fH=figure
    subplot(2,2,1)
    minNumPts=30;
    %minNumPts=40;
    delPhaseSmoothed(dataCountPerBinXTP<minNumPts)=NaN;
    delSpeedSmoothed(dataCountPerBinXTS<minNumPts)=NaN;
    delCycleDursSmoothed(dataCountPerBinXTC<minNumPts)=NaN;
   omarPcolor(delDistBins,delTimeBins,-(delPhaseSmoothed'-180)/360,fH)
        %xlabel('Distance along track (cm)')
       if(normEachCell)
             xlabel('\Deltax (fraction of field)')
             ylabel('\Deltat (fraction of field)')
        else
            xlabel('\Deltax (cm)')
             ylabel('\Deltat (sec)')
        end
        colormap(jet)
        cb=colorbar
        %ylabel(cb,'Mean \Deltap precession (degrees)')
        ylabel(cb,'Mean \Deltap precession (fraction of theta cycle)')
       %ylim([0 maxDispTime])
      %caxis([0 360])
       %daspect([1 1 1])
       hold on
        %caxis([0 270])
        %caxis([0 60])
        caxis([0 1/5.5])
        caxis([0 1/4.5])
        caxis([0 1/6])
         %caxis([0.025 0.15])
          %caxis([0.07 0.125])
          %caxis([0.06 0.15])
          caxis([0.1 0.25])
         
        %caxis([50 120])
        %caxis([40 120])
      xlim([minDistXTP maxDistXTP])
      %ylim([0 0.25])
       ylim([minTimeXTP maxTimeXTP])
       hold on
       plot(xlim,ylim,'k--','LineWidth',3)
       
       title(sprintf('Average theta phase change in delta-cycle spacetime'))
        daspect([1 1 1])
             subplot(2,2,2)
          omarPcolor(delDistBins,delTimeBins,dataCountPerBinXTP',fH)
        %xlabel('Distance along track (cm)')
        if(normEachCell)
             xlabel('\Deltax (fraction of field)')
             ylabel('\Deltat (fraction of field)')
        else
            xlabel('\Deltax (cm)')
             ylabel('\Deltat (sec)')
        end
        colormap(gca,parula)
        cb=colorbar
        %ylabel(cb,'Mean \Deltap precession (degrees)')
        ylabel(cb,'Count')
       %ylim([0 maxDispTime])
      %caxis([0 360])
      %caxis([0 180])
       %caxis([0 100])
        %caxis([0 1/4])
      xlim([minDistXTP maxDistXTP])
      %ylim([0 0.25])
      ylim([minTimeXTP maxTimeXTP])
       daspect([1 1 1])
       hold on
       plot(xlim,ylim,'k--','LineWidth',3)
       
       subplot(2,2,3)
          omarPcolor(delDistBinsCycleDurs,delTimeBinsCycleDurs,1./delCycleDursSmoothed',fH)
        %xlabel('Distance along track (cm)')
        if(normEachCell)
             xlabel('\Deltax (fraction of field)')
             ylabel('\Deltat (fraction of field)')
        else
            xlabel('\Deltax (cm)')
             ylabel('\Deltat (sec)')
        end
        colormap(gca,parula)
        cb=colorbar
        caxis([6.5 8])
        ylabel(cb,'LFP theta frequency (Hz)')

      xlim([minDistXTP maxDistXTP])
      ylim([minTimeXTP maxTimeXTP])
       daspect([1 1 1])
       hold on
       plot(xlim,ylim,'k--','LineWidth',3)
       subplot(2,2,4)
          omarPcolor(delDistBinsSpeed,delTimeBinsSpeed,delSpeedSmoothed',fH)
        %xlabel('Distance along track (cm)')
        if(normEachCell)
             xlabel('\Deltax (fraction of field)')
             ylabel('\Deltat (fraction of field)')
        else
            xlabel('\Deltax (cm)')
             ylabel('\Deltat (sec)')
        end
        colormap(gca,jet)
        cb=colorbar
        ylabel(cb,'Speed (cm/s)')
        caxis([10 35])
        caxis([0 40])

      xlim([minDistXTP maxDistXTP])
      ylim([minTimeXTP maxTimeXTP])
       daspect([1 1 1])
       hold on
       plot(xlim,ylim,'k--','LineWidth',3)
       
       maxFig
       setFigFontTo(18)

   
       print('-r75',sprintf('deltaSpacetimeVsPhaseChangeSummary_%s',codeVersionID),'-dpng')
       
       %plotSpacetimePredictedPhaseCone
       %close all

end

%another property of log compression is linear relationship between time and width 
close all; clc

clearvars -except data

 %base=1.4;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load time and phase data within good field traversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setTightSubplots_Spacious
setPastelColors

plotConditionalProb=1;
%plotConditionalProb=0;

showPlots=0;

linearFitBinnedAvgs=1;
usePeakPhase=1;

numTimeBins=50;
numTimeBins=20;
numTimeBins=15;

numTimeBins=25;
numTimeBins=30;
%numTimeBins=50;
%numTimeBins=35;

timeMin=0.025;
%timeMin=0;

timeMin=0.15;
timeMin=0.05;
timeMin=0.1;
timeMin=0.01;
timeMin=0;
if(timeMin==0)
    timeMinLog=0;
else
    timeMinLog=getLogTime(timeMin);
end
numTimeAvgBins=10;
numTimeAvgBins=30;
%numTimeAvgBins=50;
%timeFieldAvgEdges=linspace(timeMin,1,numTimeAvgBins+1);
timeFieldAvgEdges=linspace(0,1,numTimeAvgBins+1);

timeFieldEdges=linspace(timeMin,1,numTimeBins+1);
logTimeFieldEdges=linspace(timeMinLog,1,numTimeBins+1);

numPhaseBins=50;
numPhaseBins=120;
numPhaseBins=60;
%numPhaseBins=30;
%numPhaseBins=90;
phaseEdges=linspace(0,360,numPhaseBins+1);
phaseErrEdges=linspace(-180,180,numPhaseBins+1);

 %bases=[1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0];
 
 %bases=1.1:0.05:2;
 
  bases=1.05:0.05:2;
  
    bases=1.05:0.05:1.7;
   
 
 totalAbsoluteErrorPerBaseDistTimeLdistLtime=NaN(length(bases),4);
 
 for b=1:length(bases)
     
         base=bases(b);
        
        

        %intermediateFileName='Feb9_2022_phaseVsTimeAndLogTimeStats.mat';
        intermediateFileName='Feb10_2022_phaseVsTimeAndLogTimeStats.mat';
        if(~exist(intermediateFileName,'file'))

            %data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
            tic
            disp('loading spike data....')
            if(~exist('data','var'))
                data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
            end
            toc

            saveFieldPhaseVsLogTimeDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3Images/fieldPhaseVsLogTimeDir';
            touchDir(saveFieldPhaseVsLogTimeDir)

            saveFigureDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/figurePanels';
            touchDir(saveFigureDir)

            %totalUsedFieldCount=data.totalUsedFieldCount;
            totalFieldCount=data.totalFieldCount;

            allFieldSpikePhases=data.allUnitSpikePhases;

            %spikeTimeFracInFieldPerField=data.allUnitSpikeTimesInField;
            %{
            allRenormedSpikeTimeFracs=[];
            for fi=1:totalFieldCount
                currFieldRenormedSpikes=data.autoRenormedSpikeTimeFracInFieldPerField{fi};

                allRenormedSpikeTimeFracs=[allRenormedSpikeTimeFracs; currFieldRenormedSpikes(:)];
            end
            spikeTimeFracInFieldPerField=allRenormedSpikeTimeFracs;
            %}

            allSpikeTimeFracs=[];
            allSpikeDistFracs=[];

            for fi=1:totalFieldCount

                currFieldSpikeDistFracs=data.distInFieldPerSpikePerField{fi};
                currFieldSpikeTimeFracs=data.spikeTimeFracInFieldPerField{fi};

                allSpikeDistFracs=[allSpikeDistFracs; currFieldSpikeDistFracs(:)];
                allSpikeTimeFracs=[allSpikeTimeFracs; currFieldSpikeTimeFracs(:)];
            end
            spikeTimeFracInFieldPerField=allSpikeTimeFracs;
            spikeDistFracInFieldPerField=allSpikeDistFracs;



                %allFieldToFieldEndTimes=1-spikeTimeFracInFieldPerField;

                allFieldToFieldEndTimes=spikeTimeFracInFieldPerField;
                allFieldToFieldEndDists=spikeDistFracInFieldPerField;


                 %base=1.5;

                timeOffset=0.125;
                 timeOffset=0.05;
                 timeOffset=0;

                allFieldToFieldEndLogTimes=1-getLogTime(1-spikeTimeFracInFieldPerField,base);
                allFieldToFieldEndLogDists=1-getLogTime(1-spikeDistFracInFieldPerField,base);

                %{
                allFieldToFieldEndLogTimes(allFieldToFieldEndLogTimes<=0)=NaN;
                 allFieldToFieldEndLogTimes(allFieldToFieldEndLogTimes>=1)=NaN;


                 allFieldToFieldEndLogDists(allFieldToFieldEndLogDists<=0)=NaN;
                 allFieldToFieldEndLogDists(allFieldToFieldEndLogDists>=1)=NaN;
                %}

                %allFieldToFieldEndLogTimes=getLogTime(spikeTimeFracInFieldPerField,base);



                if(showPlots)
                    fH=figure;
                    %title(sprintf('Phase vs log(time to end of field), all fields'))
                    %xlabel('log(time to end of field)')
                    %ylabel('Theta phase (deg)')
                end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %phase vs behavior stats
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                     if(linearFitBinnedAvgs)
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         %phase vs log(dist)
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         if(usePeakPhase)
                             [logDistBinCenters,logDistPhaseAvgValues]=...
                                 getBinnedCircPeakVals(allFieldToFieldEndLogDists(:),allFieldSpikePhases(:),timeFieldAvgEdges);
                         else
                              [logDistBinCenters,logDistPhaseAvgValues,lowerLimPerBin,upperLimPerBin]=...
                                     getBinnedCircAverages(allFieldToFieldEndLogDists(:),allFieldSpikePhases(:),timeFieldAvgEdges);
                         end


                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         %phase vs (dist)
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           if(usePeakPhase)
                             [distBinCenters,distPhaseAvgValues]=...
                                 getBinnedCircPeakVals(allFieldToFieldEndDists(:),allFieldSpikePhases(:),timeFieldAvgEdges);
                           else
                              [distBinCenters,distPhaseAvgValues,lowerLimPerBin,upperLimPerBin]=...
                                 getBinnedCircAverages(allFieldToFieldEndDists(:),allFieldSpikePhases(:),timeFieldAvgEdges);
                           end

                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         %phase vs log(time)
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         if(usePeakPhase)
                             [logTimeBinCenters,logTimePhaseAvgValues]=...
                                 getBinnedCircPeakVals(allFieldToFieldEndLogTimes(:),allFieldSpikePhases(:),timeFieldAvgEdges);
                         else
                           [logTimeBinCenters,logTimePhaseAvgValues,lowerLimPerBin,upperLimPerBin]=...
                              getBinnedCircAverages(allFieldToFieldEndLogTimes(:),allFieldSpikePhases(:),timeFieldAvgEdges);
                         end

                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         %phase vs time
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       if(usePeakPhase)
                         [timeBinCenters,timeBinPhaseAvgValues]=...
                             getBinnedCircPeakVals(allFieldToFieldEndTimes(:),allFieldSpikePhases(:),timeFieldAvgEdges);
                       else
                          [timeBinCenters,timeBinPhaseAvgValues,lowerLimPerBin,upperLimPerBin]=...
                                   getBinnedCircAverages(allFieldToFieldEndTimes(:),allFieldSpikePhases(:),timeFieldAvgEdges);
                       end


                         %[m,b,R,p]=getLinearFit(x,y)
                         linLinFit=0;
                         if(linLinFit)
                             [mLogDist,bLogDist,RLogDist,~]=getLinearFit(logDistBinCenters,logDistPhaseAvgValues);
                            [mBehavDist,bBehavDist,RBehavDist,~]=getLinearFit(distBinCenters,distPhaseAvgValues);
                            [mLogTime,bLogTime,RLogTime,~]=getLinearFit(logTimeBinCenters,logTimePhaseAvgValues);
                            [mBehavTime,bBehavTime,RBehavTime,~]=getLinearFit(timeBinCenters,timeBinPhaseAvgValues);
                         else
                            [RLogDist,pLogDist,sLogDist,bLogDist ]=kempter_lincirc(logDistBinCenters,2*pi/360*logDistPhaseAvgValues);
                            [RBehavDist,pBehavDist,sBehavDist,bBehavDist ]=kempter_lincirc(distBinCenters,2*pi/360*distPhaseAvgValues);
                            [RLogTime,pLogTime,sLogTime,bLogTime ]=kempter_lincirc(logTimeBinCenters,2*pi/360*logTimePhaseAvgValues);
                            [RBehavTime,pBehavTime,sBehavTime,bBehavTime ]=kempter_lincirc(timeBinCenters,2*pi/360*timeBinPhaseAvgValues);

                            mLogDist=sLogDist*360; %cycles/field --> degrees/field
                            bLogDist=mod(bLogDist*360/(2*pi),360); %offset, convert radians to degrees

                            mBehavDist=sBehavDist*360; %cycles/field --> degrees/field
                            bBehavDist=mod(bBehavDist*360/(2*pi),360); %offset, convert radians to degrees

                            mLogTime=sLogTime*360; %cycles/field --> degrees/field
                            bLogTime=mod(bLogTime*360/(2*pi),360); %offset, convert radians to degrees

                            mBehavTime=sBehavTime*360; %cycles/field --> degrees/field
                            bBehavTime=mod(bBehavTime*360/(2*pi),360); %offset, convert radians to degrees

                         end
                         %%
                         if(showPlots)
                             figure; 
                             yLims=[30 330];
                             yLims=[0 360];
                             markerStyle='.';
                             markerSize=30;
                             markerColors=(timeWarpPaperColorsColdToHot);

                             subplot(1,4,3)
                             plot(logDistBinCenters,logDistPhaseAvgValues,markerStyle,'MarkerSize',markerSize,'Color',markerColors(3,:));hold on; plot([0 1],[bLogDist bLogDist+mLogDist],'k--','LineWidth',3)
                             xlabel('log(Dist) (norm)')
                             ylabel('Peak phase (deg)')
                             title('Log distance')
                             axis square
                             box off
                             ylim(yLims)
                              subplot(1,4,1)
                             plot(distBinCenters,distPhaseAvgValues,markerStyle,'MarkerSize',markerSize,'Color',markerColors(1,:));hold on; plot([0 1],[bBehavDist bBehavDist+mBehavDist]+360,'k--','LineWidth',3)
                              xlabel('dist in field (norm)')
                             ylabel('Peak phase (deg)')
                             title('Linear distance')
                             axis square
                             box off

                             ylim(yLims)
                             subplot(1,4,4)
                             plot(logTimeBinCenters,logTimePhaseAvgValues,markerStyle,'MarkerSize',markerSize,'Color',markerColors(4,:));hold on; plot([0 1],[bLogTime bLogTime+mLogTime],'k--','LineWidth',3)
                             ylim(yLims)
                             xlabel('log(time) (norm)')
                             ylabel('Peak phase (deg)')
                             title('Log time')
                             axis square
                             box off
                              subplot(1,4,2)
                             plot(timeBinCenters,timeBinPhaseAvgValues,markerStyle,'MarkerSize',markerSize,'Color',markerColors(2,:));hold on; plot([0 1],[bBehavTime bBehavTime+mBehavTime],'k--','LineWidth',3)
                            xlabel('time in field (norm)')
                             ylabel('Peak phase (deg)')
                             title('Linear time')
                             ylim(yLims)
                             axis square
                             box off

                             maxFig
                             setFigFontTo(16)

                             saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,'linearVsLog_DistVsTime_PhaseLinearFit.png'))

                             disp('')
                         end

                     end


                if(showPlots)
                    cProbLims=[0 14e-4];
                    fTimeH=figure; 
                end 
                    plotPhaseVsTimeAndLogTime

               if(showPlots)
                    setFigFontTo(18)
                    maxFig

                   saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('logTimeVsThetaPhase_AllFields.png')))

                     fDistH=figure; 

               end
                    plotPhaseVsDistAndLogDist


                if(showPlots)

                    %set(gca,'xscale','log')
                    setFigFontTo(18)
                    maxFig

                    saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('logDistVsThetaPhase_AllFields.png')))
                    close all

                end

               %plot mean error with confidence interval vs field progression
               %%

             numErrBins=10;
             timeFieldAvgErrEdges=linspace(0,1,numErrBins);
             usePeakError=1;

               if(~usePeakError)
                   [timeBinCenters,timeBinPhaseErrAvgValues,lowerLimErrPerTimeBin,upperLimErrPerTimeBin]=...
                                   getBinnedCircAverages(allFieldToFieldEndTimes(:),behavTimeErrors(:),timeFieldAvgErrEdges);

                    [logTimeBinCenters,logTimePhaseErrAvgValues,lowerLimErrPerLogTimeBin,upperLimErrPerLogTimeBin]=...
                                  getBinnedCircAverages(allFieldToFieldEndLogTimes(:),logTimeErrors(:),timeFieldAvgErrEdges);

                   [distBinCenters,distPhaseErrAvgValues,lowerLimErrPerDistBin,upperLimErrPerDistBin]=...
                                     getBinnedCircAverages(allFieldToFieldEndDists(:),behavDistErrors(:),timeFieldAvgErrEdges);

                     [logDistBinCenters,logDistPhaseErrAvgValues,lowerLimErrPerLogDistBin,upperLimErrPerLogDistBin]=...
                             getBinnedCircAverages(allFieldToFieldEndLogDists(:),logDistErrors(:),timeFieldAvgErrEdges);
                           ylabelStr='Mean error (deg)';

               else

                 [timeBinCenters,timeBinPhaseErrAvgValues]=...
                           getBinnedCircPeakVals(allFieldToFieldEndTimes(:),behavTimeErrors(:),timeFieldAvgErrEdges);

                [logTimeBinCenters,logTimePhaseErrAvgValues]=...
                              getBinnedCircPeakVals(allFieldToFieldEndLogTimes(:),logTimeErrors(:),timeFieldAvgErrEdges);

               [distBinCenters,distPhaseErrAvgValues]=...
                                 getBinnedCircPeakVals(allFieldToFieldEndDists(:),behavDistErrors(:),timeFieldAvgErrEdges);

                 [logDistBinCenters,logDistPhaseErrAvgValues]=...
                         getBinnedCircPeakVals(allFieldToFieldEndLogDists(:),logDistErrors(:),timeFieldAvgErrEdges);
                    ylabelStr='Peak phase error (deg)';
               end
               
               currBaseTotalDistAbsErr=sum(abs(distPhaseErrAvgValues));
               currBaseTotalTimeAbsErr=sum(abs(timeBinPhaseErrAvgValues));
               currBaseTotalLogDistAbsErr=sum(abs(logDistPhaseErrAvgValues));
               currBaseTotalLogTimeAbsErr=sum(abs(logTimePhaseErrAvgValues));

                     figure; 

                         yLims=[-180 180];
                         yLims=[-50 50];
                         %yLims=[-100 100];
                         markerStyle='.';
                         markerSize=30;
                         markerColors=(timeWarpPaperColorsColdToHot);

                         subplot(1,4,3)
                         plot(logDistBinCenters,logDistPhaseErrAvgValues,markerStyle,'MarkerSize',markerSize,'Color',markerColors(3,:));
                         hold on
                         hline(0,'k--',3)
                         xlabel(sprintf('log_{%.2f}(distance) (norm)',base))
                          ylabel(ylabelStr)
                         title(sprintf('log_{%.2f}(distance)',base))
                         axis square
                         box off
                         ylim(yLims)
                          subplot(1,4,1)
                         plot(distBinCenters,distPhaseErrAvgValues,markerStyle,'MarkerSize',markerSize,'Color',markerColors(1,:));
                          xlabel('dist in field (norm)')
                         ylabel(ylabelStr)
                         title('Linear distance')
                         axis square
                         box off
                         hold on
                         hline(0,'k--',3)

                         ylim(yLims)
                         subplot(1,4,4)
                         plot(logTimeBinCenters,logTimePhaseErrAvgValues,markerStyle,'MarkerSize',markerSize,'Color',markerColors(4,:));
                         ylim(yLims)
                         xlabel(sprintf('log_{%.2f}(time) (norm)',base))
                          ylabel(ylabelStr)
                         title(sprintf('log_{%.2f}(time)',base))
                         axis square
                         box off
                         hold on
                         hline(0,'k--',3)

                          subplot(1,4,2)
                         plot(timeBinCenters,timeBinPhaseErrAvgValues,markerStyle,'MarkerSize',markerSize,'Color',markerColors(2,:));
                        xlabel('time in field (norm)')
                          ylabel(ylabelStr)
                         title('Linear time')
                         ylim(yLims)
                         axis square
                         box off

                         hold on 
                         hline(0,'k--',3)

                         maxFig
                         setFigFontTo(16)

                         totalAbsoluteErrorPerBaseDistTimeLdistLtime(b,1)=currBaseTotalDistAbsErr;
                         totalAbsoluteErrorPerBaseDistTimeLdistLtime(b,2)=currBaseTotalTimeAbsErr;
                         totalAbsoluteErrorPerBaseDistTimeLdistLtime(b,3)=currBaseTotalLogDistAbsErr;
                         totalAbsoluteErrorPerBaseDistTimeLdistLtime(b,4)=currBaseTotalLogTimeAbsErr;
                         
                         saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('linearVsLogBase%.2f_DistVsTime_PhaseErrLinearFit.png',base)))




            %%

            %save('phaseVsTimeAndLogTimeStats.mat')

        else
            load(intermediateFileName)
            %load('phaseVsTimeAndLogTimeStats.mat')
        end
    
 end
 %%
 figure
 for m=1:4
     plot(bases,totalAbsoluteErrorPerBaseDistTimeLdistLtime(:,m),'Color',markerColors(m,:),'LineWidth',3)
     hold on
 end
 axis tight
 ylim([0 max(totalAbsoluteErrorPerBaseDistTimeLdistLtime(:))])
 xlabel('Logarithm base')
 ylabel('Total error (deg)')
 legend('linear dist','linear time','log dist','log time','Location','northwest')
 box off
 legend boxoff
 
%saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('linearVsLogBaseAbsErr_DistVsTime_PhaseErrLinearFit.png')))



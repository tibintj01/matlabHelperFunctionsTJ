%another property of log compression is linear relationship between time and width 
close all; clc

clearvars -except data
setPastelColors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load time and phase data within good field traversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setTightSubplots_Spacious

plotConditionalProb=1;
%plotConditionalProb=0;


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

    showPlots=1;
    
        %allFieldToFieldEndTimes=1-spikeTimeFracInFieldPerField;
        
        allFieldToFieldEndTimes=spikeTimeFracInFieldPerField;
        allFieldToFieldEndDists=spikeDistFracInFieldPerField;
        
        %{
        %domain restriction also in logTime, to make fair comparison!
        allFieldToFieldEndTimes(allFieldToFieldEndTimes<=0)=NaN;
         allFieldToFieldEndTimes(allFieldToFieldEndTimes>=1)=NaN;
         
         
         allFieldToFieldEndDists(allFieldToFieldEndDists<=0)=NaN;
         allFieldToFieldEndDists(allFieldToFieldEndDists>=1)=NaN;
        %}

         %base=1.475;
         %base=sqrt(2);
          %base=1.33333333;
         base=1.4;
         %base=1.375;
         
         %base=1.45;
        %allFieldToFieldEndLogTimes=getLogTime(1-spikeTimeFracInFieldPerField-0.15,base);
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
                 
             else
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %phase vs log(dist)
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                notNaNIdxes=~isnan(allFieldToFieldEndLogDists) & ~isnan(allFieldSpikePhases) & timeMinLog<=allFieldToFieldEndLogDists;
                [RLogDist,pLogDist,sLogDist,bLogDist ] = kempter_lincirc(allFieldToFieldEndLogDists(notNaNIdxes),2*pi/360*allFieldSpikePhases(notNaNIdxes));
                mLogDist=sLogDist*360; %cycles/field --> degrees/field
                bLogDist=mod(bLogDist*360/(2*pi),360); %offset, convert radians to degrees
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %phase vs (dist)
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 notNaNIdxes=~isnan(allFieldToFieldEndDists) & ~isnan(allFieldSpikePhases) & timeMin<=allFieldToFieldEndDists;
                [RBehavDist,pDist,sDist,bDist ] = kempter_lincirc(allFieldToFieldEndDists(notNaNIdxes),2*pi/360*allFieldSpikePhases(notNaNIdxes));
                mBehavDist=sDist*360; %cycles/field --> degrees/field
                bBehavDist=mod(bDist*360/(2*pi),360); %offset, convert radians to degrees
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %phase vs log(time)
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                notNaNIdxes=~isnan(allFieldToFieldEndLogTimes) & ~isnan(allFieldSpikePhases) & timeMinLog<=allFieldToFieldEndLogTimes;
               [RLogTime,p,s,b ] = kempter_lincirc(allFieldToFieldEndLogTimes(notNaNIdxes),2*pi/360*allFieldSpikePhases(notNaNIdxes));
                mLogTime=s*360; %cycles/field --> degrees/field
                bLogTime=mod(b*360/(2*pi),360); %offset, convert radians to degrees

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %phase vs time
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 notNaNIdxes=~isnan(allFieldToFieldEndTimes) & ~isnan(allFieldSpikePhases) & timeMin<=allFieldToFieldEndTimes;
                [RBehavTime,p,s,b ] = kempter_lincirc(allFieldToFieldEndTimes(notNaNIdxes),2*pi/360*allFieldSpikePhases(notNaNIdxes));
                mBehavTime=s*360; %cycles/field --> degrees/field
                bBehavTime=mod(b*360/(2*pi),360); %offset, convert radians to degrees
             end
            
          
        if(showPlots)
            cProbLims=[0 14e-4];
            fTimeH=figure; 
            plotPhaseVsTimeAndLogTime
            setFigFontTo(18)
            maxFig
            
           saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('logTimeVsThetaPhase_AllFields.png')))

             fDistH=figure; 
            plotPhaseVsDistAndLogDist
            
            %set(gca,'xscale','log')
            setFigFontTo(18)
            maxFig
            
            saveas(gcf,fullfile(saveFieldPhaseVsLogTimeDir,sprintf('logDistVsThetaPhase_AllFields.png')))


            close all

        end

       


    

    %%

    %save('phaseVsTimeAndLogTimeStats.mat')

else
    load(intermediateFileName)
    %load('phaseVsTimeAndLogTimeStats.mat')
end
    
%%
%numFieldsUsed
setTightSubplots_Spacious

numRhoBins=50;
rhoEdges=linspace(-1,1,numRhoBins+1);
rhoBins=edgesToBins(rhoEdges);

numErrorBins=50;
numErrorBins=100;
minError=-90;
maxError=90;

minError=-80;
maxError=80;
errorEdges=linspace(minError,maxError,numErrorBins+1);

errorBins=edgesToBins(errorEdges);

figure
subplot(2,1,1)
Nlog=histcounts(allFieldLogTimePhaseRs,rhoEdges);
Nlin=histcounts(allFieldTimePhaseRs,rhoEdges);
pLinRho=plot(rhoBins,smooth(Nlin),'k-','LineWidth',4)
hold on
pLogRho=plot(rhoBins,smooth(Nlog),'r-','LineWidth',4)

title(sprintf('circR distributions, %d fields',numFieldsUsed ))
xlim([-1 1])

legend([pLinRho pLogRho],{'phase-time linear model','phase-log(time) linear model'},'Location','northwest')

xlabel('Circular correlation coefficient')
ylabel('Field count')

box off
subplot(2,1,2)
NlogErr=histcounts(allFieldLogTimeMeanCircErrors,errorEdges);
NlinErr=histcounts(allFieldTimeMeanCircErrors,errorEdges);
pLinErr=plot(errorBins,smooth(NlinErr),'k-','LineWidth',4)
hold on
pLogErr=plot(errorBins,smooth(NlogErr),'r-','LineWidth',4)
xlabel('Mean error (deg)')
ylabel('Field count')
legend([pLinErr pLogErr],{'phase-time linear model','phase-log(time) linear model'},'Location','northwest')
box off

title(sprintf('error distributions, %d fields',numFieldsUsed ))

setFigFontTo(18)
maxFigHalfWidth


saveas(gcf,'linearModelLogTimeFitStatsDistr.png')
%xlim([-1 1])

%{
histogram(allFieldLogTimePhaseSlopes,30)
title(sprintf('circSlope distribution (deg/field), %d fields',numFieldsUsed ))

xlim([-360 720])
%}






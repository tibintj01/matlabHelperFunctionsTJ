close all; clear all; clc
%dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';

figurePanelDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/figurePanels';
touchDir(figurePanelDir)

allDataByField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');


processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';
unitDataDir=fullfile(processedDataDir,'unitSpikeInfoStructs');
dataDir=unitDataDir;

if(~exist(fullfile(processedDataDir,'allPhaseJumpData.mat'),'file'))

    filePaths=getFilePathsRegex(dataDir,'*mat');

    minNumPts=5;
    minNumThetaCyclesInField=3;
    %minNumThetaCyclesInField=10;

    minSeqSpeed=0.05; %m/s
    %minSeqSpeed=0.4; %m/s
    %minSeqSpeed=0.3; %m/s
    %minSeqSpeed=0.6; %m/s

    %maxSeqSpeed=0.8; %m/s
    %minSeqSpeed=0.45; %m/s
    allPhaseJumpRs=[];
    allPhaseJumpSlopes=[];

    masterTimesInField=[];
    masterStartPhaseBinCenters=[];
    masterTimeBinCenters=[];
    masterPhaseGains=[];
    masterPhaseGainsMovAvg=[];
    startPhaseLims=[50 310];
    startPhaseLims=[0 360];
    %startPhaseLims=[90 270];
    numTimeBins=10;
    timeFieldEdges=linspace(0,1,numTimeBins+1);

    timeWindowFields=0.1;
    movingAvgStepSize=0.01;
    movingAvgStepSize=0.025;
    
    
    movingAvgStartPhaseStepSize=5;
    movingAvgStartPhaseWindHalf=12.5;
    movingAvgStartPhaseWindHalf=30;
    %movingAvgStartPhaseWindHalf=25;
    movingAvgStartPhaseWindHalf=12.5;

    %movingAvgWindHalf=0.05;
    %movingAvgWindHalf=0.1;
    %movingAvgWindHalf=0.075;
    movingAvgWindHalf=0.025;
    
    %movingAvgWindHalf=0.1;
    movAvgTimeBinCenters=movingAvgStepSize:(movingAvgStepSize):(1-movingAvgStepSize);
    movAvgStartPhaseBinCenters=movingAvgStartPhaseStepSize:(movingAvgStartPhaseStepSize):(360-movingAvgStartPhaseStepSize);

    useLinCircR=1;
    useLinCircR=0;

    useInvSlope=1;
    useInvSlope=0;

    numPlots=0;

    showPlots=0;

    %%
    fM=figure(11)

    justSummary=1;

    totalFieldCount=allDataByField.totalFieldCount;

    numFieldsUsedTotal=0;
    timeJitterSigma=0.05;
    timeJitterSigma=0.025;
timeJitterSigma=0;

    for fi=1:totalFieldCount
           fi
        currFilePath=allDataByField.unitInfoPathPerField{fi};

        currFileName=getFileNameFromPath(currFilePath);
        fileBaseName=currFileName(1:(end-4));

        data=load(currFilePath);


                currPhaseJumpR=NaN;
                currPhaseJumpSlope=NaN;



                 currTimesInField=allDataByField.timeInFieldPerCyclePerField{fi};
                 currPhaseGains=allDataByField.phaseGainsPerCyclePerField{fi}./allDataByField.timeGainsPerCyclePerField{fi};
                 
                 
                 
                 %outOfBoundsIdxes=isnan(currTimesInField) | currTimesInField>1 | currPhaseGains>180 | currPhaseGains<-180;
                 outOfBoundsIdxes=isnan(currTimesInField) | currTimesInField>1;


                 
                 currTimesInField(outOfBoundsIdxes)=NaN;
                 currPhaseGains(outOfBoundsIdxes)=NaN;
                 
                 %{
                 thetaJumpFracDiff=max(diff((unique(currTimesInField(~isnan(currTimesInField))))));
                 if(isempty(thetaJumpFracDiff))
                     continue
                 end
                 
                 
                 %currPhaseGains=currPhaseGains/thetaJumpFracDiff;
                 
                 approxNumThetaCyclesInField=1/thetaJumpFracDiff;
                 
                 if(approxNumThetaCyclesInField<minNumThetaCyclesInField) %so that cycle represents gain
                     continue
                 end
                 
                 %}
                
                 currTimesInField=currTimesInField+normrnd(0,timeJitterSigma,size(currTimesInField));

                 if(length(currPhaseGains)<minNumPts)
                     continue
                 end

                 notNaNIdxes=~isnan(currTimesInField) & ~isnan(currPhaseGains);
                [phaseGainMovAvgValues] = getBinnedMovingAverages(currTimesInField(notNaNIdxes),currPhaseGains(notNaNIdxes),movingAvgWindHalf,movAvgTimeBinCenters);

                 %if(~justSummary)
                     %figure
                     %plot(currTimesInField,currPhaseGains,'k.','MarkerSize',30)

                     if(showPlots)
                         
                         %{
                          if(useInvSlope)
                              plot(currTimesInField,1./currPhaseGains,'k.','MarkerSize',10)
                          else
                             plot(currTimesInField,currPhaseGains,'k.','MarkerSize',10)
                          end
                         %}
                          
                          plot(movAvgTimeBinCenters,phaseGainMovAvgValues,'k.','MarkerSize',10)

                          xlim([0 1])

                          %plot(currStartPhasesInField,1./currPhaseGains,'k.','MarkerSize',30)
                         hold on
                     end
                % end


                if(useInvSlope)
                    currPhaseGains=1./currPhaseGains;
                end
                 %rMat=corrcoef([phaseJumpVsTimeVsDistInField(corrIdxes,2),phaseJumpVsTimeVsDistInField(corrIdxes,3)]);
                 
                 
                 notNaNIdxesAvg=~isnan(phaseGainMovAvgValues) & ~isnan(movAvgTimeBinCenters);
                 %rMat=corrcoef([currTimesInField(notNaNIdxes),currPhaseGains(notNaNIdxes)]);
                 %rMat=corrcoef([movAvgTimeBinCenters(notNaNIdxesAvg),phaseGainMovAvgValues(notNaNIdxesAvg)]);

                 [m,b,R]=getLinearFit(movAvgTimeBinCenters(notNaNIdxesAvg),phaseGainMovAvgValues(notNaNIdxesAvg));


                 currPhaseJumpSlope=m;
                 currPhaseJumpR=R;
                 %{
                 if(length(rMat)>1)
                        currPhaseJumpR=rMat(1,2);
                 %linCoeffsJump=polyfit(phaseJumpVsTimeVsDistInField(corrIdxes,2),phaseJumpVsTimeVsDistInField(corrIdxes,3),1);
                              linCoeffsJump=polyfit(currTimesInField(notNaNIdxes),currPhaseGains(notNaNIdxes),1);
                            currPhaseJumpSlope=linCoeffsJump(1);
                 else
                     currPhaseJumpR=NaN;
                     currPhaseJumpSlope=NaN;
                 end

                 %}

                 if(useLinCircR)
                     currPhaseJumpR
                     linVar=currPhaseGains;
                     circVar=currStartPhasesInField*pi/180;
                     [ currPhaseJumpR,p,s,b ] = kempter_lincirc( linVar,circVar );
                     currPhaseJumpSlope=s*360;
                 end




                 %if(useInvSlope)
                 %        phaseGainMovAvgValues=1./phaseGainMovAvgValues;
                 %    end



                 %[phaseGainMovAvgValues] = getBinnedMovingAverages(currTimesInField,currPhaseGains,movingAvgStartPhaseWindHalf,movAvgStartPhaseBinCenters)

                 %{
                 [sortedTimesInField,sortedTimesI]=sort(currTimesInField(:));
                 timeSortedPhaseGains=currPhaseGains(sortedTimesI);


                 %[timeBinCenters,phaseGainAverages] = getBinnedAverages(currTimesInField,currPhaseGains,timeFieldEdges)
                 %plot(timeBinCenters,phaseGainAverages,'r-','LineWidth',5)

                 approxFieldsPerIdx=max(0.01,nanmedian(diff(sortedTimesInField)));
                 timeWindIdx=round(timeWindowFields/approxFieldsPerIdx);
                 phaseGainSmoothed=movmean(timeSortedPhaseGains,timeWindIdx);
                 plot(sortedTimesInField,phaseGainSmoothed,'r-','LineWidth',5)
                 hold on
                 plot(sortedTimesInField,phaseGainSmoothed,'ko','MarkerSize',5)
                 %}


                 if(~justSummary)
                     %plot(movAvgTimeBinCenters,phaseGainMovAvgValues,'r-','LineWidth',5)
                      %plot(movAvgStartPhaseBinCenters,phaseGainMovAvgValues,'r-','LineWidth',5)
                 else
                           figure(fM)
                 end

                     hold on
                     %plot(movAvgTimeBinCenters,phaseGainMovAvgValues,'k.','MarkerSize',5)




                          %plot(movAvgTimeBinCenters,phaseGainMovAvgValues,'k.','MarkerSize',8)
                      numPlots=numPlots+1;

              masterTimeBinCenters=[masterTimeBinCenters;  movAvgTimeBinCenters(:)];
               %masterStartPhaseBinCenters=[masterStartPhaseBinCenters; movAvgStartPhaseBinCenters(:)];
               masterPhaseGainsMovAvg=[masterPhaseGainsMovAvg; phaseGainMovAvgValues(:)];

                 masterTimesInField=[masterTimesInField;  currTimesInField(:)];
               %masterStartPhaseBinCenters=[masterStartPhaseBinCenters; movAvgStartPhaseBinCenters(:)];
               masterPhaseGains=[masterPhaseGains; currPhaseGains(:)];

                allPhaseJumpRs=[allPhaseJumpRs currPhaseJumpR];
                allPhaseJumpSlopes=[allPhaseJumpSlopes currPhaseJumpSlope];

                numFieldsUsedTotal=numFieldsUsedTotal+1;

                %{
                if(abs(currPhaseJumpR)>0.99)
                    disp('')
                end
                %}

                if(~justSummary)
                        %title({removeUnderscores(fileBaseName),'cycle-to-cycle phase jump vs time in field'})
                      title({removeUnderscores(fileBaseName),'cycle-to-cycle phase jump vs initial phase'})

                    xlim([0 1])
                    ylim([-180 0])
                    xlabel('time in field (frac)')
                    ylabel('phase jump to next cycle (degrees)')
                    setFigFontTo(18)
                    %saveas(gcf,sprintf('%s_PhaseJumpVsTime_%s.tif',fileBaseName,currFieldDirection))
                    daspect([1 1 1])
                    saveas(gcf,sprintf('%s_PhaseJumpVsStartPhase.tif',fileBaseName))
                    close all
                end  

    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %mean of individual cell mean phase gains vs time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [masterPhaseGainMovAvgValues,masterPhaseGainMovSEMValues] = getBinnedMovingAverages(masterTimeBinCenters,masterPhaseGains,movingAvgWindHalf,movAvgTimeBinCenters);
    %[masterPhaseGainMovAvgValues,masterPhaseGainMovSEMValues] = getBinnedMovingAverages(masterStartPhaseBinCenters,masterPhaseGains,movingAvgStartPhaseWindHalf,movAvgStartPhaseBinCenters);

    save(fullfile(processedDataDir,'allPhaseJumpData.mat'))
else
    load(fullfile(processedDataDir,'allPhaseJumpData.mat'))
    
end
    
try
    figure(fM)
catch
    figure
end

plot(masterTimesInField,masterPhaseGains,'k.')
hold on

%plot(movAvgTimeBinCenters,masterPhaseGainMovAvgValues,'r-','LineWidth',5)
%hold on
             %plot(movAvgTimeBinCenters,masterPhaseGainMovAvgValues,'ko','MarkerSize',5)
             
             %sH=shadedErrorBar(movAvgTimeBinCenters,masterPhaseGainMovAvgValues,masterPhaseGainMovSEMValues,'lineprops',{'r-','linewidth',3})
              %sH=shadedErrorBar(movAvgStartPhaseBinCenters,masterPhaseGainMovAvgValues,masterPhaseGainMovSEMValues,'lineprops',{'r-','linewidth',3})
              startRegressionPhase=30;
              endRegressionPhase=330;
              startRegressionPhase=50;
              endRegressionPhase=310;
              
              %startRegressionPhase=10;
              %endRegressionPhase=350;
              nonNanIdx= ~isnan(masterPhaseGains) & ~isnan(masterTimesInField);
                linCoeffsJump=polyfit(masterTimesInField(nonNanIdx),masterPhaseGains(nonNanIdx),1);
                linModelStartTimes=[0 1];
                linModelPhaseGains=linCoeffsJump(1)*linModelStartTimes+linCoeffsJump(2);
                
                rMat=corrcoef([masterTimesInField(nonNanIdx),masterPhaseGains(nonNanIdx)]);
                
                Rval=rMat(1,2);
                %{
                    Rval = 0.5295
                
                %}
                
              linRegH=plot(linModelStartTimes,linModelPhaseGains,'r-','LineWidth',5)
              hold on
             
             xlabel('Time in field (frac)')
             %xlabel('Initial phase (deg)')
             ylabel('Cycle-to-cycle phase gain (deg)')
             

             setFigFontTo(28)
                         
            hold on
            pH=plot(0,0,'k.')
            %{
                 %ylim([-110 -40])
                  ylim([-150 -20])
                 
                 
                 ylim([-180 0])
                 %xlim(startPhaseLims)
                 xlim([0 330])
                     xlim([0 350])
                     xlim([0 360])
                     xlim([0 720])
                      xlim([-360 720])
            %}
                 
                  %lgd=legend([sH.mainLine pH], {sprintf('mean%sSEM (n=%d)',char(177),length(filePaths)),'individual cells'})
                  %lgd=legend([linRegH pH], {sprintf('linear regression (n=%d)',length(filePaths)),'individual cells'})
                  %lgd=legend([linRegH pH], {sprintf('linear regression (n=%d cells)',length(filePaths))},'Location','Southeast')
                  lgd=legend([linRegH pH], {sprintf('linear regression (n=%d fields)',numFieldsUsedTotal)},'Location','Southeast')
                  
                  %lgd=legend('individual cells')
           lgd.FontSize=12;
                     % daspect([1 .5 1])
                      %daspect([1 .25 1])
                      
                      legend boxoff
                 maxFigHalfWidth
                 
                    xlim([0 1])
             %daspect([1 25 1])
             %uberTitle(sprintf('Population average phase gain vs time in field (min speed: %.2f m/s)',minSeqSpeed))
             uberTitle(sprintf('HC3 population temporal phase gain vs starting phase (min speed: %.2f m/s)',minSeqSpeed))

          
             
             saveas(gcf,fullfile(figurePanelDir,'fig3C_GainVsTime_HC3.tif'))
             print(fullfile(figurePanelDir,'fig3C_GainVsTime_HC3'),'-depsc')
             

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distribution of phase gain slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
rEdges=-0.6:0.05:0.6;
rEdges=-0.6:0.025:0.6;
rEdges=-0.8:0.025:0.8;

rEdges=-1:0.025:1;
%rEdges=-0.6:0.00625:0.6;
%rEdges=-0.6:0.1:0.6;
slopeEdges=-75:10:75;
slopeEdges=-75:8:75;
slopeEdges=-75:6:75;
slopeEdges=-75:3:75;

slopeEdges=-100:3:100;

slopeEdges=-1000:3:1000;
%{
slopeEdges=-.3:.02:.3;
slopeEdges=-1:.025:1;
%}
%slopeEdges=-2.5:.025:2.5;
subplot(2,1,1)
[N,edges]=histcounts(allPhaseJumpRs,rEdges);
Nsmoothed=smooth(N/sum(N));
plot(edgesToBins(rEdges),Nsmoothed,'k-','LineWidth',5)
hold on
plot([0 0], ylim,'k--','LineWidth',5)
%legend(sprintf('n=%d cell-direction pairs',length(allPhaseJumpRs(~isnan(allPhaseJumpRs)))))
%legend(sprintf('n=%d cells',length(filePaths)))
legend(sprintf('n=%d cells',numFieldsUsedTotal))

legend boxoff
xlabel('cycle-to-cycle phase gain vs phase, R')
ylabel('probability')
box off

%{
[p,h,stats] = signrank(allPhaseJumpRs)

    p = 5.2912e-10
    zval: 6.2102
    signedrank: 1417

    nanmean(allPhaseJumpRs)=0.3480
    getSEMacrossRows(allPhaseJumpRs')= 0.0264

%}

subplot(2,1,2)
[N,edges]=histcounts(allPhaseJumpSlopes,slopeEdges);
Nsmoothed=smooth(N/sum(N));


plot(edgesToBins(slopeEdges),Nsmoothed,'k-','LineWidth',5)
hold on

yLims=ylim;
plot([0 0], yLims,'k--','LineWidth',5)
ylim(yLims)
box off
%legend(sprintf('n=%d cell-direction pairs',length(allPhaseJumpSlopes(~isnan(allPhaseJumpSlopes)))))
legend(sprintf('n=%d fields',numFieldsUsedTotal))

legend boxoff
%xlabel('phase gain vs time slope [(deg shift/cycle)/field]')
xlabel('cycle-to-cycle phase gain vs time, slope')
ylabel('probability')

setFigFontTo(18)

%{
[p,h,stats] = signrank(allPhaseJumpSlopes)

    p = 8.2838e-10
    zval: 6.1394
    signedrank: 1409



nanmean(allPhaseJumpSlopes)= 0.1949
getSEMacrossRows(allPhaseJumpSlopes')= 0.0241


%}

uberTitle(sprintf('HC3 phase gain vs time correlations (min speed: %.2f m/s)',minSeqSpeed))

saveas(gcf,fullfile(figurePanelDir,'fig3D_GainVsTimeDistr_HC3.tif'))
             print(fullfile(figurePanelDir,'fig3D_GainVsTimeDistr_HC3'),'-depsc')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%average phase gain across population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





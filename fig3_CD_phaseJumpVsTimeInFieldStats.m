close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
filePaths=getFilePathsRegex(dataDir,'*mat');


minSeqSpeed=0.05; %m/s
%minSeqSpeed=0.4; %m/s
%minSeqSpeed=0.3; %m/s
%minSeqSpeed=0.6; %m/s

%maxSeqSpeed=0.8; %m/s
%minSeqSpeed=0.45; %m/s
allPhaseJumpRs=[];
allPhaseJumpSlopes=[];

masterTimeBinCenters=[];
masterStartPhaseBinCenters=[];
masterPhaseGains=[];
startPhaseLims=[50 310];
startPhaseLims=[0 360];
%startPhaseLims=[90 270];
numTimeBins=10;
timeFieldEdges=linspace(0,1,numTimeBins+1);

timeWindowFields=0.1;
movingAvgStepSize=0.01;
movingAvgStartPhaseStepSize=5;
movingAvgStartPhaseWindHalf=12.5;
movingAvgStartPhaseWindHalf=30;
%movingAvgStartPhaseWindHalf=25;
movingAvgStartPhaseWindHalf=12.5;

movingAvgWindHalf=0.05;
%movingAvgWindHalf=0.075;
%movingAvgWindHalf=0.025;
%movingAvgWindHalf=0.1;
movAvgTimeBinCenters=movingAvgStepSize:(movingAvgStepSize):(1-movingAvgStepSize);
movAvgStartPhaseBinCenters=movingAvgStartPhaseStepSize:(movingAvgStartPhaseStepSize):(360-movingAvgStartPhaseStepSize);

useLinCircR=1;
useLinCircR=0;

useInvSlope=1;
useInvSlope=0;

numPlots=0;

%%
fM=figure(11)

justSummary=1;

for fi=1:length(filePaths)
       fi
    currFilePath=filePaths{fi};
    currFileName=getFileNameFromPath(currFilePath);
    fileBaseName=currFileName(1:(end-4));
    
    data=load(currFilePath);
    
    for di=1:2
        
        if(di==2)
            currFieldDirection='leftward';     
        elseif(di==1)
            currFieldDirection='rightward'; 
        end
            currPhaseJumpR=NaN;
            currPhaseJumpSlope=NaN;
            if(di==1 && isfield(data,'rPhaseJumpVsTimeRight'))
                %currPhaseJumpR=data.rPhaseJumpVsTimeRight;
                %currPhaseJumpSlope=data.jumpVsTimeSlopeRight;
                phaseJumpVsTimeVsDistInField=data.phaseJumpVsTimeVsDistInFieldRightward;
            elseif(di==2 && isfield(data,'rPhaseJumpVsTimeLeft'))
                %currPhaseJumpR=data.rPhaseJumpVsTimeLeft;
                %currPhaseJumpSlope=data.jumpVsTimeSlopeLeft;
                phaseJumpVsTimeVsDistInField=data.phaseJumpVsTimeVsDistInFieldLeftward;
            else
                continue
            end

           if(isnan(phaseJumpVsTimeVsDistInField))
               continue
           end
            
             corrIdxes=~isnan(phaseJumpVsTimeVsDistInField(:,2)) & ~isnan(phaseJumpVsTimeVsDistInField(:,3)) & phaseJumpVsTimeVsDistInField(:,2)<=1; %only within typical field duration
             corrIdxes=corrIdxes & phaseJumpVsTimeVsDistInField(:,6)>=minSeqSpeed;
             %& phaseJumpVsTimeVsDistInField(:,5)<=maxSeqSpeed;
             
             currTimesInField=phaseJumpVsTimeVsDistInField(corrIdxes,2);
             
            currStartPhasesInField=phaseJumpVsTimeVsDistInField(corrIdxes,7);
             currPhaseGains=phaseJumpVsTimeVsDistInField(corrIdxes,3);
             
             
             if(~justSummary)
                 figure
                 %plot(currTimesInField,currPhaseGains,'k.','MarkerSize',30)
                 plot(currStartPhasesInField,currPhaseGains,'k.','MarkerSize',30)

                  %plot(currStartPhasesInField,1./currPhaseGains,'k.','MarkerSize',30)
                 hold on
             end
             
             
             %rMat=corrcoef([phaseJumpVsTimeVsDistInField(corrIdxes,2),phaseJumpVsTimeVsDistInField(corrIdxes,3)]);
                     rMat=corrcoef([currStartPhasesInField(:),currPhaseGains(:)]);
                     
                     

             if(length(rMat)>1)
                    currPhaseJumpR=rMat(1,2);
             %linCoeffsJump=polyfit(phaseJumpVsTimeVsDistInField(corrIdxes,2),phaseJumpVsTimeVsDistInField(corrIdxes,3),1);
                          linCoeffsJump=polyfit(currStartPhasesInField,currPhaseGains,1);
                        currPhaseJumpSlope=linCoeffsJump(1);
             else
                 currPhaseJumpR=NaN;
                 currPhaseJumpSlope=NaN;
             end
        
       
             if(useLinCircR)
                 currPhaseJumpR
                 linVar=currPhaseGains;
                 circVar=currStartPhasesInField*pi/180;
                 [ currPhaseJumpR,p,s,b ] = kempter_lincirc( linVar,circVar );
                 currPhaseJumpSlope=s*360;
             end
       
           

             %[phaseGainMovAvgValues] = getBinnedMovingAverages(currTimesInField,currPhaseGains,movingAvgWindHalf,movAvgTimeBinCenters)
             [phaseGainMovAvgValues] = getBinnedMovingAverages(currStartPhasesInField,currPhaseGains,movingAvgStartPhaseWindHalf,movAvgStartPhaseBinCenters)

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
                 
                 if(useInvSlope)
                     phaseGainMovAvgValueInv=1./phaseGainMovAvgValues;
                     plot(movAvgStartPhaseBinCenters-360,phaseGainMovAvgValueInv,'k.','MarkerSize',8)
                   hold on
                    plot(movAvgStartPhaseBinCenters+360-movingAvgStartPhaseStepSize-360,phaseGainMovAvgValueInv,'k.','MarkerSize',8)
                   plot(movAvgStartPhaseBinCenters+720-movingAvgStartPhaseStepSize*2-360,phaseGainMovAvgValueInv,'k.','MarkerSize',8)
                   %ylim([-0.15 0])
                   ylim([-0.05 0])
                 else
                   plot(movAvgStartPhaseBinCenters-360,phaseGainMovAvgValues,'k.','MarkerSize',8)
                   hold on
                    plot(movAvgStartPhaseBinCenters+360-movingAvgStartPhaseStepSize-360,phaseGainMovAvgValues,'k.','MarkerSize',8)
                   plot(movAvgStartPhaseBinCenters+720-movingAvgStartPhaseStepSize*2-360,phaseGainMovAvgValues,'k.','MarkerSize',8)
                   numPlots=numPlots+1;
                   %xlim([-40 400])
                   %subplot(11,11,numPlots)

                 end

           
           masterTimeBinCenters=[masterTimeBinCenters;  movAvgTimeBinCenters(:)];
           masterStartPhaseBinCenters=[masterStartPhaseBinCenters; movAvgStartPhaseBinCenters(:)];
           masterPhaseGains=[masterPhaseGains; phaseGainMovAvgValues(:)];
           
            allPhaseJumpRs=[allPhaseJumpRs currPhaseJumpR];
            allPhaseJumpSlopes=[allPhaseJumpSlopes currPhaseJumpSlope];
            
            if(~justSummary)
                    %title({removeUnderscores(fileBaseName),'cycle-to-cycle phase jump vs time in field'})
                  title({removeUnderscores(fileBaseName),'cycle-to-cycle phase jump vs initial phase'})

                xlim([0 360])
                ylim([-180 0])
                xlabel('initial phase (degrees)')
                ylabel('phase jump to next cycle (degrees)')
                setFigFontTo(18)
                %saveas(gcf,sprintf('%s_PhaseJumpVsTime_%s.tif',fileBaseName,currFieldDirection))
                daspect([1 1 1])
                saveas(gcf,sprintf('%s_PhaseJumpVsStartPhase_%s.tif',fileBaseName,currFieldDirection))
                close all
            end  
    end 
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mean of individual cell mean phase gains vs time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[masterPhaseGainMovAvgValues,masterPhaseGainMovSEMValues] = getBinnedMovingAverages(masterTimeBinCenters,masterPhaseGains,movingAvgWindHalf,movAvgTimeBinCenters);
[masterPhaseGainMovAvgValues,masterPhaseGainMovSEMValues] = getBinnedMovingAverages(masterStartPhaseBinCenters,masterPhaseGains,movingAvgStartPhaseWindHalf,movAvgStartPhaseBinCenters);

try
    figure(fM)
catch
    figure
end
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
              nonNanIdx=~isnan(masterStartPhaseBinCenters) & ~isnan(masterPhaseGains) & masterStartPhaseBinCenters>=startRegressionPhase &  masterStartPhaseBinCenters<=endRegressionPhase;
                linCoeffsJump=polyfit(masterStartPhaseBinCenters(nonNanIdx),masterPhaseGains(nonNanIdx),1);
                linModelStartPhases=[0 360];
                linModelPhaseGains=linCoeffsJump(1)*linModelStartPhases+linCoeffsJump(2);
                
                rMat=corrcoef([masterStartPhaseBinCenters(nonNanIdx),masterPhaseGains(nonNanIdx)]);
                
                Rval=rMat(1,2);
                %{
                    Rval = 0.5295
                
                %}
                
              linRegH=plot(linModelStartPhases,linModelPhaseGains,'r-','LineWidth',5)
              hold on
             
             %xlabel('Time in field (frac)')
             xlabel('Initial phase (deg)')
             ylabel('Cycle-to-cycle phase gain (deg)')
             

             setFigFontTo(28)
                         
            hold on
            pH=plot(0,0,'k.')
                 %ylim([-110 -40])
                  ylim([-150 -20])
                 
                 
                 ylim([-180 0])
                 %xlim(startPhaseLims)
                 xlim([0 330])
                     xlim([0 350])
                     xlim([0 360])
                     xlim([0 720])
                      xlim([-360 720])
                 
                  %lgd=legend([sH.mainLine pH], {sprintf('mean%sSEM (n=%d)',char(177),length(filePaths)),'individual cells'})
                  %lgd=legend([linRegH pH], {sprintf('linear regression (n=%d)',length(filePaths)),'individual cells'})
                  lgd=legend([linRegH pH], {sprintf('linear regression (n=%d cells)',length(filePaths))},'Location','Southeast')
                  %lgd=legend('individual cells')
           lgd.FontSize=12;
                      daspect([1 .5 1])
                      daspect([1 .25 1])
                      
                      legend boxoff
                 maxFigHalfWidth
             %daspect([1 25 1])
             %uberTitle(sprintf('Population average phase gain vs time in field (min speed: %.2f m/s)',minSeqSpeed))
             uberTitle(sprintf('Population average temporal phase gain vs starting phase (min speed: %.2f m/s)',minSeqSpeed))


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distribution of phase gain slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
rEdges=-0.6:0.05:0.6;
rEdges=-0.6:0.025:0.6;
rEdges=-0.8:0.025:0.8;
%rEdges=-0.6:0.00625:0.6;
%rEdges=-0.6:0.1:0.6;
slopeEdges=-75:10:75;
slopeEdges=-75:8:75;
slopeEdges=-75:6:75;
slopeEdges=-75:3:75;
slopeEdges=-.3:.02:.3;
slopeEdges=-1:.025:1;
%slopeEdges=-2.5:.025:2.5;
subplot(2,1,1)
[N,edges]=histcounts(allPhaseJumpRs,rEdges);
Nsmoothed=smooth(N/sum(N));
plot(edgesToBins(rEdges),Nsmoothed,'k-','LineWidth',5)
hold on
plot([0 0], ylim,'k--','LineWidth',5)
%legend(sprintf('n=%d cell-direction pairs',length(allPhaseJumpRs(~isnan(allPhaseJumpRs)))))
legend(sprintf('n=%d cells',length(filePaths)))
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
legend(sprintf('n=%d cells',length(filePaths)))

legend boxoff
%xlabel('phase gain vs time slope [(deg shift/cycle)/field]')
xlabel('cycle-to-cycle phase gain vs phase, slope')
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

uberTitle(sprintf('Phase gain vs phase correlations (min speed: %.2f m/s)',minSeqSpeed))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%average phase gain across population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





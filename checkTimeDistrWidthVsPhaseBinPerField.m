%another property of log compression is linear relationship between time and width 
close all; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load time and phase data within good field traversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exampleChoices=[359,1094,1197];

setTightSubplots_Medium
data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';

saveDataPath=fullfile(processedDataDir,'phaseDistrWidthVsTimeInFieldData.mat');

saveFieldLogCompressPhaseVarDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/fieldLogCompressPhaseVarDir';
touchDir(saveFieldLogCompressPhaseVarDir)

saveFigureDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/figurePanels';
touchDir(saveFigureDir)

%totalUsedFieldCount=data.totalUsedFieldCount;
totalFieldCount=data.totalFieldCount;

spikePhasesPerField=data.spikePhasesPerField;
unitInfoPaths=data.unitInfoPathPerField;
spikeTimeFracInFieldPerField=data.spikeTimeFracInFieldPerField;
%spikeTimeFracInFieldPerField=data.spikeTimeFracInFieldPerFieldRenormed;



useKernelEstDistr=1;
useStdDev=1;


numTimeBins=200;

maxTimeWidth=1;
maxTime=1.3;
minTime=-0.5;
maxTime=1.5;
maxTime=1.3;
maxCalcTime=1.3;
maxTime=1.2;
maxCalcTime=1.2;

maxTime=1;
maxCalcTime=1;

timeEdges=linspace(minTime,maxTime,numTimeBins+1);
%timeEdges=linspace(0,1.5,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);
timeBinWidth=median(diff(timeBinCenters));



gaussKernelWidth=0.1;
gaussKernelWidth=0.15;
gaussKernelWidth=0.1;
%gaussKernelWidth=0.05;

numPhaseBins=7;
%numPhaseBins=10;

phaseEdges=linspace(0,360,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);
phaseBinWidth=median(angdiffDeg(phaseBinCenters));

phaseBinColors=jet(numPhaseBins);


minNumSpikesForField=300;
minNumSpikesForField=200;
minNumSpikesForField=100;
%minNumSpikesForField=50;

minNumSpikesPerPhaseBin=10;
%minNumSpikesPerPhaseBin=5;
%minNumSpikesPerPhaseBin=3;
%minNumSpikesPerPhaseBin=20;

numFieldsUsed=0;

showPlots=0;
%showPlots=1;

allFieldTimeWidthSlopes=[];
allFieldTimeWidthRs=[];
allFieldTimeWidthPs=[];
allFieldTimeBinCenters=[];
allFieldDistrWidths=[];
allFieldLinearModelErrs=[];

allFieldTimeMeanRange=[];


allFieldStatIDs=[];
allFieldTimeWidthAvgs=[];


useDistrEstimation=0;

for fi=1:totalFieldCount
    
%for fi=exampleChoices
    currFieldSpikePhases=spikePhasesPerField{fi};
    unitInfoStruct=load(unitInfoPaths{fi});
    sessionName=unitInfoStruct.unitInfo.sessionName;
    
    %currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi};
    
     if(length(currFieldSpikePhases)<=minNumSpikesForField)
        continue
    else
          
        numFieldsUsed=numFieldsUsed+1
        %continue
    end
    
    %currFieldInFieldTimes=getRealTimeBoundNormedSpikeFracs(spikePhasesPerField{fi},spikeTimeFracInFieldPerField{fi});
    
       %currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi}*1.15;
    
     currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi}/0.8; %actual time bound

     if(size(currFieldInFieldTimes,2)>1)
            currFieldInFieldTimes=currFieldInFieldTimes(:,1);
        end
     %currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi}; %actual time bound
     
     currFieldInFieldTimes(currFieldInFieldTimes>maxCalcTime)=NaN;
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %shift until zero phase has circmean 360
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     [bestShiftDeg,shiftedPhases]=getBestShiftDegZeroPhase(currFieldSpikePhases,currFieldInFieldTimes);
          %[bestShiftDeg,shiftedPhases]=getBestShiftDeg1D(currFieldSpikePhases,currFieldInFieldTimes,0.25);

     
          %{
     figure; 
     subplot(2,1,1); plot(currFieldInFieldTimes,currFieldSpikePhases,'ko')
     subplot(2,1,2); plot(currFieldInFieldTimes,shiftedPhases,'ko')
     currFieldSpikePhases=shiftedPhases;
          %}
     
    
     %if(showPlots)
     %     fH=figure;
     %     [pxySmooth,pxy]=getJointDistr(currFieldInFieldTimes,currFieldSpikePhases,timeEdges,phaseEdges,fH);
     %else
     if(useDistrEstimation)
          [pxySmooth,pxy]=getJointDistr(currFieldInFieldTimes,currFieldSpikePhases,timeEdges,phaseEdges);
     end
     %end

   
    
    maxProbTimePerPhaseBin=NaN(numPhaseBins,1);
    distrWidth=NaN(numPhaseBins,1);
    correspondingTimeBinCenters=NaN(numPhaseBins,1);
    

    for phi=1:(numPhaseBins)
        m=NaN;
        R=NaN;
        b=NaN;
        
        %get kappa for current bin
         currPhaseBinStart=phaseBinCenters(phi)-phaseBinWidth/2;
            currPhaseBinEnd=phaseBinCenters(phi)+phaseBinWidth/2;
            
            inBinIdxes=currFieldSpikePhases>=currPhaseBinStart & currFieldSpikePhases<currPhaseBinEnd;
            
           inBinPhases=currFieldSpikePhases(inBinIdxes);
           inBinTimes=currFieldInFieldTimes(inBinIdxes);
           
           if(length(inBinPhases)<minNumSpikesPerPhaseBin || sum(~isnan(inBinTimes))==0 )
               continue
           end
           
            notNaNIdxesPhase=~isnan(inBinPhases);
            
       
            if(useDistrEstimation)
            
            if(useKernelEstDistr)
                
                 kernelDistrTimeObj = fitdist(inBinTimes(:),'Kernel','BandWidth',gaussKernelWidth); %gaussian kernel estimation, 1/10th of field
                  %phaseAxis=linspace(0,360,361);
                 %phaseAxis=linspace(0,360,numPhaseBins+1);
                 currPdist=pdf(kernelDistrTimeObj,timeBinCenters);
                 %currPdist=circSmooth(currPdist,30);
                
               
            else
                assert(false) %stop code, this condition not implemented
                %currPdist=circSmooth(pxySmooth(phi,:),3);
            end
            %currPdist=circSmooth(pxy(ti,:),5);
            %currPdist=circSmooth(pxy(ti,:),5);

            currPdist=currPdist/sum(currPdist);


              [~,maxID]=max(currPdist);
        %    currPdist=circshift(currPdist,-maxID);

           
             
           %[currDistrWidth,risingIntersectPhase,risingIntersectProb,descendingIntersectPhase,descendingIntersectProb]...
           %    =getFullWidthFracMaxCircDistr(phaseBinCenters,currPdist,0.5);
           
           %{
           [currDistrWidth,risingIntersectTime,risingIntersectProb,descendingIntersectTime,descendingIntersectProb]...
               =getFullWidthHalfMaxLinearDistr(timeBinCenters,currPdist);
           
           correspondingTimeBinCenter=(risingIntersectTime+descendingIntersectTime)/2;
           
           if(currDistrWidth>maxTime)
               continue
           end
           %}
            end
           
            noOutliersPerPhaseBin=1;
            minNumOutlierDetect=6;
            outlierAlpha=0.05;
            sameSize=1;
           if(useStdDev)
               if(noOutliersPerPhaseBin)
                   if(sum(~isnan(inBinTimes))>=minNumOutlierDetect)
                       inBinTimesNoOutliers=deleteoutliers(inBinTimes,outlierAlpha,sameSize);
                       %{
                       figure;
                       pO=plotLineHist(inBinTimes,linspace(0,1,50),'k');
                       hold on
                       pNO=plotLineHist(inBinTimesNoOutliers,linspace(0,1,50),'r');
                       legend([pO pNO],{'outliers','nooutliers'});
                       %}
                       inBinTimes=inBinTimesNoOutliers;
                   end
               end
               
               currDistrWidth=nanstd(inBinTimes);
               correspondingTimeBinCenter=nanmean(inBinTimes);
           end
           
        if(showPlots)
            if(useDistrEstimation)
                subplot(1,3,1)
            else
                subplot(1,2,1)
            end
            
            %plot(currFieldInFieldTimes,currFieldSpikePhases,'k.','MarkerSize',25)
            plot(currFieldInFieldTimes,currFieldSpikePhases,'k.','MarkerSize',20)
            
            hold on
            
            
            
            ylim([0 360])

            

           if(useDistrEstimation)
                subplot(1,3,2)
            
       
            plot(timeBinCenters,currPdist,'-','Color',phaseBinColors(phi,:),'LineWidth',3)

            hold on

            if(~(useStdDev))
            plot(risingIntersectTime,risingIntersectProb,'.','Color',phaseBinColors(phi,:),'MarkerSize',30)
            plot(descendingIntersectTime,descendingIntersectProb,'.','Color',phaseBinColors(phi,:),'MarkerSize',30)

            plot([risingIntersectTime descendingIntersectTime], [risingIntersectProb descendingIntersectProb],'Color',phaseBinColors(phi,:),'LineWidth',3)
            end

            hold on
           end
           
           if(useDistrEstimation)
              subplot(1,3,3)
           else
               subplot(1,2,2)
           end
            
          
                
                plot(correspondingTimeBinCenter,currDistrWidth,'.','Color',phaseBinColors(phi,:),'MarkerSize',80)
   
     
            hold on
           
            
            
            
         end
             
                 distrWidth(phi)=currDistrWidth;
                 correspondingTimeBinCenters(phi)=correspondingTimeBinCenter;
             
            %maxProbTimePerPhaseBin(phi)=timeBinCenters(maxID);
    end
    
    notNaNIdxes=~isnan(distrWidth) & distrWidth<=maxTimeWidth & correspondingTimeBinCenters<=maxTimeWidth;
        %notNaNIdxes=~isnan(distrWidth); %& distrWidth<=maxTimeWidth & correspondingTimeBinCenters<=maxTimeWidth;

        correspondingTimeBinCentersNonNan=correspondingTimeBinCenters(notNaNIdxes);
        distrWidthNonNan=distrWidth(notNaNIdxes);
        
        currFieldAvgTimeWidth=mean(distrWidthNonNan);
        %currFieldAvgPhaseWidth=0;
       %currFieldAvgPhaseWidth=nanmin(distrWidthNonNan);
       
       currFieldTimeMeanRange=range(correspondingTimeBinCentersNonNan);
        
        [m,b,R,p]=getLinearFit(correspondingTimeBinCentersNonNan,distrWidthNonNan);
        
        %{
        if(exp(-m)>1.5 || exp(-m)<1.3)
            showPlots=0;
        else
            showPlots=1;
        end
        
         if(exp(-m)>1.5 || exp(-m)<1.3)
             close all
                continue
        end
        %}
        
            y0=b;
            yf=m+b;
            
        %currFieldAvgPhaseWidth=b;
        
        predictedVals=interp1([minTime maxTime],[y0 yf],correspondingTimeBinCentersNonNan);
        actualVals=distrWidthNonNan;

        signedLinErrorsRealMinusPredicted=(actualVals(:)-predictedVals(:));
        
       

        %timeBinCenters
        
        allFieldTimeBinCenters=[allFieldTimeBinCenters;correspondingTimeBinCentersNonNan(:)];
        allFieldDistrWidths=[allFieldDistrWidths;distrWidthNonNan(:)];      
        allFieldLinearModelErrs=[allFieldLinearModelErrs;signedLinErrorsRealMinusPredicted(:)];
        
        allFieldStatIDs=[allFieldStatIDs;repelem(numFieldsUsed,length(distrWidthNonNan),1)];

    if(showPlots)
        
            if(useDistrEstimation)
                subplot(1,3,1)
            else
                subplot(1,2,1)
            end
         hold on
         for bi=1:numPhaseBins
            errorbar(correspondingTimeBinCenters(bi),phaseBinCenters(bi),distrWidth(bi),'horizontal','.','MarkerSize',60,'Color',phaseBinColors(bi,:),'LineWidth',4,'CapSize',20)
         end
    xlabel('Time in field (frac)')
    ylabel('Theta phase (deg)')
    xlim([0 maxCalcTime])
    xlim([-0.1 maxCalcTime+0.1])
    ylim([0 360])
        title({sessionName,sprintf('Theta phase vs Time in traversal, field %d (n=%d spikes)',fi,length(currFieldSpikePhases))})

        colormap(jet)
     cb1=colorbar;
     caxis([0 360])
     ylabel(cb1,'Theta phase bin center (deg)')
      box off
     
         if(useDistrEstimation)
            subplot(1,3,2)
            colormap(gca,jet)
            cb=colorbar
            ylabel(cb,'Theta phase (deg)')
            xlabel('Time in field (frac)')

            ylabel('Probability')
            %xlim([0 360])
            %xlim([minTime maxTime])
            xlim([-0.1 maxCalcTime+0.1])
            %xlim([0 maxTime])
            caxis([0 360])

            title(sprintf('Time in field distribution vs theta phase bin'))
         end


    %daspect([1 360 1])
    %caxis([0 prctile(pxySmooth(:),97.5)])
    %maxFig
    if(useDistrEstimation)
        subplot(1,3,3)
    else
        subplot(1,2,2)
    end
    
    
       
       title({sprintf('Distribution width vs time, R=%.2f',R),sprintf('slope=%.2f, estimated log base: %.2f',m,exp(-m))})
    
    
        hold on
    %x0=minTime;
    x0=0;
    y0=m*x0+b;
    xf=maxTime;
    yf=m*xf+b;
    plot([x0 xf], [y0 yf],'k--','LineWidth',5)
    
     xlabel('Mean time for phase bin (field frac)')

            ylabel('Std time for phase bin (field frac)')
    xlim([0 maxTime])
    %ylim([0 1])
     %ylim([0 0.5])
      ylim([0 0.4])
     colormap(jet)
     cb2=colorbar;
     caxis([0 360])
     ylabel(cb2,'Theta phase bin center (deg)')

     axis square
     box off
     setFigFontTo(18)
    maxFig
    saveas(gcf,fullfile(saveFieldLogCompressPhaseVarDir,sprintf('timeVsTimeWidthForThetaPhase_Session_%s_Base%.2f_Field%d.png',sessionName,exp(-m),fi)))
    %saveEPS(fullfile(saveFieldLogCompressPhaseVarDir,sprintf('timeVsTimeWidthForThetaPhase_Field%d.eps',fi)))
    close all

    
    end
    
    allFieldTimeWidthSlopes=[allFieldTimeWidthSlopes; m];
    allFieldTimeWidthRs=[allFieldTimeWidthRs; R];
    allFieldTimeWidthPs=[allFieldTimeWidthPs; p];
    allFieldTimeWidthAvgs=[allFieldTimeWidthAvgs; currFieldAvgTimeWidth];
    
    allFieldTimeMeanRange=[allFieldTimeMeanRange; currFieldTimeMeanRange];

end

%%
%numFieldsUsed

figure

plotLineHist((allFieldTimeWidthRs),linspace(-1,1,51))
xlabel('Mean vs std. dev. time R')
ylabel('Probability')
axis tight
axis square
setFigFontTo(16)
%%
figure; 
subplot(1,2,1)


slopeMean=nanmean(allFieldTimeWidthSlopes);
rThresh=0.3;
highRIdxes=(abs(allFieldTimeWidthRs)>rThresh);
rThreshSlopeMean=nanmean(allFieldTimeWidthSlopes(highRIdxes));
%getSEMacrossRows(allFieldTimeWidthSlopes(highRIdxes))
plotLineHist((-allFieldTimeWidthSlopes(highRIdxes)),linspace(0,1,51))



%histogram(exp(-allFieldTimeWidthSlopes),linspace(-1,1,31))
%histogram(exp(-allFieldTimeWidthSlopes(highRIdxes)),linspace(0,3,51))
numHighRFields=sum(highRIdxes);
xlabel('Mean vs std time estimated log base')
ylabel('Probability')

title({sprintf('Mean vs std time slope in equal width phase bins (|R|>%.2f, n=%d fields)',rThresh,numHighRFields),...
    sprintf('mean estimated log base = %.3f (e^{-m}) ',exp(-rThreshSlopeMean))})

subplot(1,2,2)
rThresh=0.25;
rThresh=0.3;
highRIdxes=(abs(allFieldTimeWidthRs)>rThresh);
rThreshSlopeMean=nanmean(allFieldTimeWidthSlopes(highRIdxes));

%fraction of fields within range
baseRangeMin=1.1;
baseRangeMax=1.7;

allFieldBases=exp(-allFieldTimeWidthSlopes(highRIdxes));
%nanmean(exp(-allFieldTimeWidthSlopes(highRIdxes)))
%>>1.4321
%getSEMacrossRows(exp(-allFieldTimeWidthSlopes(highRIdxes)))
%>> 0.0105

sum(allFieldBases>=baseRangeMin & allFieldBases<=baseRangeMax)

length(allFieldBases)

%histogram(exp(-allFieldTimeWidthSlopes),linspace(-1,1,31))
%histogram(exp(-allFieldTimeWidthSlopes(highRIdxes)),linspace(0,3,51))
plotLineHist(exp(-allFieldTimeWidthSlopes(highRIdxes)),linspace(0,3,51))
axis tight
axis square
numHighRFields=sum(highRIdxes);
xlabel('Mean vs std time estimated log base')
ylabel('Probability')

title({sprintf('Mean vs std time slope in equal width phase bins (|R|>%.2f, n=%d fields)',rThresh,numHighRFields),...
    sprintf('mean estimated log base = %.3f (e^{-m})',exp(-rThreshSlopeMean))})

setFigFontTo(16)
maxFig

saveas(gcf,'estimatedLogBaseDistrMeanVsStdTime.png')

save(saveDataPath,'allFieldStatIDs', 'allFieldTimeWidthAvgs',...
    'allFieldTimeWidthSlopes','allFieldTimeWidthRs','allFieldTimeBinCenters','allFieldDistrWidths',...
    'minNumSpikesPerPhaseBin', 'minNumSpikesForField','numFieldsUsed','timeBinCenters','allFieldTimeMeanRange','allFieldTimeWidthPs')






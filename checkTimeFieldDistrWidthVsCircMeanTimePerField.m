%another property of log compression is linear relationship between time and width 
close all; clear all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load time and phase data within good field traversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setTightSubplots_Spacious
data=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat')

processedDataDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/hc3ProcessedData';

saveDataPath=fullfile(processedDataDir,'phaseDistrWidthVsTimeInFieldData.mat')

saveFieldLogCompressPhaseVarDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/fieldLogCompressPhaseVarDir';
touchDir(saveFieldLogCompressPhaseVarDir)

saveFigureDir='/Users/tibinjohn/thetaSeq/code/hc3ProcessingCode/figurePanels';
touchDir(saveFigureDir)

%totalUsedFieldCount=data.totalUsedFieldCount;
totalFieldCount=data.totalFieldCount;

spikePhasesPerField=data.spikePhasesPerField;
spikeTimeFracInFieldPerField=data.spikeTimeFracInFieldPerField;

useKappaVarMeasure=1;
useKappaVarMeasure=0;

peakNorm=1;


useKernelEstDistr=0;
useKernelEstDistr=1;
gaussTimeKernelWidth=0.1;


estimateMeanSTDfromKernelDistr=0;
%estimateMeanSTDfromKernelDistr=1;

computeDistrWidth=0;

maxStartingTime=1/3; %cell must code for time in first third of field to ascertain whole pattern across field

takeSqrt=1;
takeSqrt=0;

baselineWidth=25;
baselineWidth=60;
%baselineWidth=75;

numTimeBins=10;
numTimeBins=8;
numTimeBins=7;
numTimeBins=10;
numTimeBins=15;
numTimeBins=10;
%numTimeBins=12;
%numTimeBins=15;
numTimeBins=15;
numTimeBins=20;

if(useKernelEstDistr)
    numTimeBins=100;
else
    numTimeBins=50;
end

%startTimeEdges=-1;
%endTimeEdges=2;

startTimeEdges=-0.5;
endTimeEdges=1.5;

%timeEdges=linspace(0,1,numTimeBins+1);
timeEdges=linspace(startTimeEdges,endTimeEdges,numTimeBins+1);
%timeEdges=linspace(0,1.5,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);
timeBinWidth=median(diff(timeBinCenters));

%minTimeDisp=-0.25;
%maxTimeDisp=1.25;

timeBinCentersDisp=(startTimeEdges+timeBinWidth/2):timeBinWidth:(endTimeEdges-timeBinWidth/2);


%numPhaseBins=24;
numPhaseBins=30;
numPhaseBins=60;
numPhaseBins=24;
numPhaseBins=45;
numPhaseBins=360;
numPhaseBins=30;
numPhaseBins=24;
%numPhaseBins=15;
numPhaseBins=8;
numPhaseBins=12;
numPhaseBins=24;
numPhaseBins=18;
numPhaseBins=12;
numPhaseBins=7;
numPhaseBins=12;
numPhaseBins=7;
numPhaseBins=18;
numPhaseBins=7;
%numPhaseBins=8;
%numPhaseBins=12;
%numPhaseBins=5;
%numPhaseBins=8;
%numPhaseBins=12;
%numPhaseBins=24;
%numPhaseBins=30;
%phaseEdges=linspace(0,360,numPhaseBins+1);

%phaseBounds=[40 320];

phaseBounds=[0 360];
phaseEdges=linspace(phaseBounds(1),phaseBounds(2),numPhaseBins+1);

phaseBinCenters=edgesToBins(phaseEdges);

phaseBinWidth=median(diff(phaseBinCenters));

phaseBinColors=jet(numPhaseBins);

minNumSpikesForField=600;
minNumSpikesForField=700;
minNumSpikesForField=300;

%minNumSpikesPerPhaseBin=10;
%minNumSpikesPerPhaseBin=6;
minNumSpikesPerPhaseBin=15;
minNumSpikesPerPhaseBin=30;
%minNumSpikesPerPhaseBin=50;


numFieldsUsed=0;

showPlots=0;
showPlots=1;

allFieldTimeWidthSlopes=[];
allFieldTimeWidthRs=[];
allFieldTimeDistrMeans=[];
allFieldTimeDistrWidths=[];
allFieldLinearModelErrs=[];
allFieldTimeWidthOffsets=[];

allFieldTimeMeansByPhaseBins=[];
allFieldTimeValsByPhaseBins=[];

allUsedFieldsTimes=[];
allUsedFieldsPhases=[];

for fi=1:totalFieldCount
    
    currFieldSpikePhases=spikePhasesPerField{fi};
    currFieldInFieldTimes=1-spikeTimeFracInFieldPerField{fi}; %TIME IN FIELD VS TIME TO END OF FIELD
    
    distrWidth=NaN(numPhaseBins,1);
    distrMean=NaN(numPhaseBins,1);
    maxProbTimePerPhaseBin=NaN(numTimeBins,1);
    
    if(length(currFieldSpikePhases)<minNumSpikesForField)
        continue

    end

    if(showPlots)
        [pxySmooth,pxy]=getJointDistr(currFieldInFieldTimes,currFieldSpikePhases,timeEdges,phaseEdges); %only needed for plots not this analysis
        figure;
        title(sprintf('Phase distribution vs time in field, field %d',fi))
        xlabel('Time in field (frac)')
        ylabel('Theta phase (deg)')
    end
    
    distrVals=cell(numPhaseBins,1);
    distrMeanRepeated=cell(numPhaseBins,1);

    allInBinPhases=[];
    allInBinTimes=[];
    for phi=1:numPhaseBins
        m=NaN;
        R=NaN;
        b=NaN;
        
        %get kappa for current bin
         currPhaseBinStart=phaseBinCenters(phi)-phaseBinWidth/2;
         currPhaseBinEnd=phaseBinCenters(phi)+phaseBinWidth/2;
            
            inBinIdxes=currFieldSpikePhases>=currPhaseBinStart & currFieldSpikePhases<currPhaseBinEnd;
            
          
           inBinTimes=currFieldInFieldTimes(inBinIdxes);
           inBinPhases=currFieldSpikePhases(inBinIdxes);
           
           %notNaNIdxesTimes=~isnan(inBinTimes) & inBinTimes<=1 & inBinTimes>=0;
             %notNaNIdxesTimes=~isnan(inBinTimes)  & inBinTimes<=1 & inBinTimes>=-0.1; %end of field is more variable
            notNaNIdxesTimes=~isnan(inBinTimes)  & inBinTimes<=1.5 & inBinTimes>=-0.5; %end of field is more variable
             %notNaNIdxesTimes=~isnan(inBinTimes); %end of field is more variable

                          
           if(sum(notNaNIdxesTimes)<=minNumSpikesPerPhaseBin)
               continue
           end
           
    
           inBinTimes=inBinTimes(notNaNIdxesTimes);
            inBinPhases=inBinPhases(notNaNIdxesTimes);
            
            allInBinTimes=[allInBinTimes; inBinTimes(:)];
            allInBinPhases=[allInBinPhases; inBinPhases(:)];
            
                
  
            %currTimeBinKappa=circ_kappa(ang2rad(inBinTimes(notNaNIdxesTimes)));
            
            if(showPlots)
                if(useKernelEstDistr)
                    kernelDistrTimeObj = fitdist(inBinTimes(:),'Kernel','BandWidth',gaussTimeKernelWidth); %gaussian kernel estimation
                    %phaseAxis=linspace(0,360,361);
                    %phaseAxis=linspace(0,360,numPhaseBins+1);
                    %currPdist=pdf(kernelDistrTimeObj,timeBinCenters);
                    currPdist=pdf(kernelDistrTimeObj,timeBinCentersDisp);
                else
                    %currPdist=circSmooth(pxySmooth(phi,:),3);
                     currPdist=circSmooth(pxySmooth(:,phi),3);
                    timeBinCentersDisp=timeBinCenters;
                end
           
            
            currPdist=currPdist/sum(currPdist);


              [~,maxID]=max(currPdist);
              if(estimateMeanSTDfromKernelDistr)
               currDistrCircMean= mean(kernelDistrTimeObj);
              end
            end
        
          
         if(~estimateMeanSTDfromKernelDistr)
             currDistrCircMean=nanmean(inBinTimes);
         end
         
         
           if(showPlots)
               if(computeDistrWidth)


                   [currDistrWidth,risingIntersectTime,risingIntersectProb,descendingIntersectTime,descendingIntersectProb]...
                       =getFullWidthHalfMaxCircDistr(timeBinCenters*360,currPdist);

                   currDistrWidth=currDistrWidth/360;

                    risingIntersectTime=risingIntersectTime/360;
                   descendingIntersectTime=descendingIntersectTime/360;

               else


                   if(estimateMeanSTDfromKernelDistr)
                    currDistrWidth=std(kernelDistrTimeObj);
                   
                   end
               end
           end
           
           if(~estimateMeanSTDfromKernelDistr)
                currDistrWidth=nanstd(inBinTimes);
           end
          
           
        if(showPlots)
          subplot(1,2,1)

          if(peakNorm)
             plot(timeBinCentersDisp,currPdist/max(currPdist(:)),'-','Color',phaseBinColors(phi,:),'LineWidth',3)
          else
               plot(timeBinCentersDisp,currPdist,'-','Color',phaseBinColors(phi,:),'LineWidth',3)
          end

            hold on

            if(computeDistrWidth)
                plot(risingIntersectTime,risingIntersectProb,'.','Color',phaseBinColors(phi,:),'MarkerSize',30)
                plot(descendingIntersectTime,descendingIntersectProb,'.','Color',phaseBinColors(phi,:),'MarkerSize',30)

                plot([risingIntersectTime descendingIntersectTime], [risingIntersectProb descendingIntersectProb],'Color',phaseBinColors(phi,:),'LineWidth',3)

            end

            %circMeanPhasePerSpeedBin(ti)=circMeanDeg(


            hold on
            subplot(1,2,2)
            
           
                
           plot(currDistrCircMean,currDistrWidth,'.','Color',phaseBinColors(phi,:),'MarkerSize',80)
   
   
            hold on
           
            %ylim()
            
           
        end
          
             distrWidth(phi)=currDistrWidth;
             distrMean(phi)=currDistrCircMean;
             distrVals{phi}=inBinTimes;
             
             distrMeanRepeated{phi}=repelem(currDistrCircMean,length(inBinTimes),1);
            
            %maxProbTimePerPhaseBin(phi)=timeBinCenters(maxID);
    end
    

            
    
    
    
        %notNaNIdxes=~isnan(distrWidth) & ~isnan(distrMean) & distrWidth<=1 & distrWidth>0 & distrMean<=1 & distrMean>0 ;
        notNaNIdxes=~isnan(distrWidth) & ~isnan(distrMean) & distrMean>=0 & distrWidth>=0;

        distrMeanNonNan=distrMean(notNaNIdxes);
        distrWidthNonNan=distrWidth(notNaNIdxes);
        
        [m,b,R]=getLinearFit(distrMeanNonNan,distrWidthNonNan);
            y0=b;
            yf=m+b;
        
        predictedVals=interp1([0 1],[y0 yf],distrMeanNonNan);
        actualVals=distrWidthNonNan;

        signedLinErrorsRealMinusPredicted=(actualVals(:)-predictedVals(:));

        %timeBinCenters
        %plot([maxStartingTime maxStartingTime],ylim,'r--','LineWidth',3)
        
        %{
        if(fi==22)
            disp('')
            
        end
        %}
        
        if(min(distrMeanNonNan)>maxStartingTime)
            close all
            continue
        end
        
        allFieldTimeDistrMeans=[allFieldTimeDistrMeans;distrMeanNonNan(:)];
        allFieldTimeDistrWidths=[allFieldTimeDistrWidths;distrWidthNonNan(:)];
        
        allFieldLinearModelErrs=[allFieldLinearModelErrs;signedLinErrorsRealMinusPredicted(:)];

    if(showPlots)
        subplot(1,2,1)
        colormap(gca,jet)
        cb=colorbar
        ylabel(cb,'Theta phase (deg)')
        xlabel('Time to end of field (frac)')

        ylabel('Probability')
        %xlim([0 360])
        %xlim([0 1])
        xlim([startTimeEdges endTimeEdges])
        caxis(phaseBounds)

        title(sprintf('Time distribution vs theta phase, field %d (n=%d spikes)',fi,length(currFieldSpikePhases)))


        %daspect([1 360 1])
        %caxis([0 prctile(pxySmooth(:),97.5)])
        %maxFig
        subplot(1,2,2)

        if(computeDistrWidth)
         title({sprintf('Time distribution width vs mean, R=%.2f',R),sprintf('slope=%.2f',m)})
        else
          title({sprintf('Time distribution std vs mean, R=%.2f',R),sprintf('slope=%.2f',m)})

        end



            hold on
        x0=0;
        y0=b;
        xf=1;
        yf=m*1+b;
        plot([x0 xf], [y0 yf],'k--','LineWidth',5)

        colormap(gca,jet)
        cb=colorbar
        ylabel(cb,'Theta phase (deg)')
        caxis(phaseBounds)

        xlabel('Time distribution mean (frac)')

         if(computeDistrWidth)
            ylabel('Time distribution width (frac)')
         else
             ylabel('Time distribution std (frac)')
         end
        %ylim([0 1])
          ylim([0 0.5])
        xlim([0 1])

        %daspect([1 1 1])
     axis square

         setFigFontTo(18)
        maxFig
        saveas(gcf,fullfile(saveFieldLogCompressPhaseVarDir,sprintf('timeVarVsMeanJointDistr_Field%d.png',fi)))
        close all

    
    end
    
    allFieldTimeWidthSlopes=[allFieldTimeWidthSlopes; m];
    allFieldTimeWidthOffsets=[allFieldTimeWidthOffsets; b];
    allFieldTimeWidthRs=[allFieldTimeWidthRs; R];
    
    
    distrValsUnfurled=[];
    distrMeanRepeatedUnfurled=[];
    
    for phi=1:numPhaseBins
        distrValsUnfurled=[distrValsUnfurled; distrVals{phi}];
        distrMeanRepeatedUnfurled=[distrMeanRepeatedUnfurled; distrMeanRepeated{phi}];
    end
    
    allFieldTimeValsByPhaseBins=[allFieldTimeValsByPhaseBins; distrValsUnfurled(:)];
    allFieldTimeMeansByPhaseBins=[allFieldTimeMeansByPhaseBins; distrMeanRepeatedUnfurled(:)];

         allUsedFieldsTimes=[allUsedFieldsTimes; allInBinTimes(:)];
    allUsedFieldsPhases=[allUsedFieldsPhases; allInBinPhases(:)];
    
    
    numFieldsUsed=numFieldsUsed+1
  

end

%%
numFieldsUsed

save(saveDataPath,'allUsedFieldsTimes','allUsedFieldsPhases', 'allFieldTimeValsByPhaseBins', 'allFieldTimeMeansByPhaseBins','allFieldTimeWidthSlopes','allFieldTimeWidthOffsets','allFieldTimeWidthRs','allFieldTimeDistrMeans','allFieldTimeDistrWidths','minNumSpikesPerPhaseBin','phaseBinCenters','allFieldLinearModelErrs')






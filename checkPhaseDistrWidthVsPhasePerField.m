%another property of log compression is linear relationship between time and width 
close all; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load time and phase data within good field traversals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%spikeTimeFracInFieldPerField=data.spikeTimeFracInFieldPerField;
spikeTimeFracInFieldPerField=data.spikeTimeFracInFieldPerFieldRenormed;
useKappaVarMeasure=1;
useKappaVarMeasure=0;

useKernelEstDistr=1;

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
numTimeBins=7;

%numTimeBins=12;
numTimeBins=15;
%numTimeBins=25;
numTimeBins=10;
numTimeBins=25;

numTimeBins=15;
numTimeBins=25;
numTimeBins=20;


timeEdges=linspace(0,1,numTimeBins+1);
%timeEdges=linspace(0,1.5,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);
timeBinWidth=median(diff(timeBinCenters));

gaussKernelWidth=60;
gaussKernelWidth=30;

gaussKernelWidth=5;

gaussKernelWidth=15;
%gaussKernelWidth=10;
%gaussKernelWidth=45;
gaussKernelWidth=30;

%numPhaseBins=24;
numPhaseBins=30;
numPhaseBins=60;
numPhaseBins=24;
numPhaseBins=30;
numPhaseBins=45;
numPhaseBins=360;
phaseEdges=linspace(0,360,numPhaseBins+1);
phaseBinCenters=edgesToBins(phaseEdges);

timeBinColors=jet(numTimeBins)


minNumSpikesForField=600;
minNumSpikesForField=700;
minNumSpikesForField=300;
%minNumSpikesForField=0;

minNumSpikesPerTimeBin=5;
minNumSpikesPerTimeBin=10;
%minNumSpikesPerTimeBin=20;
%minNumSpikesPerTimeBin=30;
%minNumSpikesPerTimeBin=30;
%minNumSpikesPerTimeBin=50;
%minNumSpikesPerTimeBin=30;
%minNumSpikesForField=3000;

numFieldsUsed=0;

showPlots=0;
%showPlots=1;

allFieldTimeWidthSlopes=[];
allFieldTimeWidthRs=[];
allFieldTimeBinCenters=[];
allFieldDistrWidths=[];
allFieldLinearModelErrs=[];
allFieldQuadModelErrs=[];

allFieldTimeWidthQuadRs=[];

allFieldStatIDs=[];
allFieldTimeWidthAvgs=[];

for fi=1:totalFieldCount
    currFieldSpikePhases=spikePhasesPerField{fi};
    
    currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi};
    
     if(length(currFieldSpikePhases)<=minNumSpikesForField)
        continue
    else
          
        numFieldsUsed=numFieldsUsed+1
        %continue
    end
    
    %currFieldInFieldTimes=getRealTimeBoundNormedSpikeFracs(spikePhasesPerField{fi},spikeTimeFracInFieldPerField{fi});
    
       %currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi}*1.15;
    
     currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi}/0.8; %actual time bound
     %currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi}; %actual time bound
    

    [pxySmooth,pxy]=getJointDistr(currFieldInFieldTimes,currFieldSpikePhases,timeEdges,phaseEdges);

    if(showPlots)
        figure;
        title(sprintf('Phase distribution vs time in field, field %d',fi))
        xlabel('Time in field (frac)')
        ylabel('Theta phase (deg)')
    end
    
    maxProbPhasePerTimeBin=NaN(numTimeBins,1);
    distrWidth=NaN(numTimeBins,1);
    

    for ti=1:numTimeBins
        m=NaN;
        R=NaN;
        b=NaN;
        
        %get kappa for current bin
         currTimeBinStart=timeBinCenters(ti)-timeBinWidth/2;
            currTimeBinEnd=timeBinCenters(ti)+timeBinWidth/2;
            
            inBinIdxes=currFieldInFieldTimes>=currTimeBinStart & currFieldInFieldTimes<currTimeBinEnd;
            
           inBinPhases=currFieldSpikePhases(inBinIdxes);
           
           if(length(inBinPhases)<minNumSpikesPerTimeBin )
               continue
           end
           
    
                
    
            notNaNIdxesPhase=~isnan(inBinPhases);
            
            currTimeBinKappa=circ_kappa(ang2rad(inBinPhases(notNaNIdxesPhase)));
            
            if(useKernelEstDistr)
                %{
                 kernelDistrPhaseObj = fitdist(inBinPhases(:),'Kernel','BandWidth',gaussKernelWidth); %gaussian kernel estimation, 30 degree
                  %phaseAxis=linspace(0,360,361);
                 %phaseAxis=linspace(0,360,numPhaseBins+1);
                 currPdist=pdf(kernelDistrPhaseObj,phaseBinCenters);
                 currPdist=circSmooth(currPdist,30);
                %}
                
                 vfObservations=inBinPhases(:)*2*pi/360;
                 vfPDFSamples=phaseBinCenters(:)*2*pi/360;
                 vfDomain=[0 2*pi];
                 fSigma=gaussKernelWidth*2*pi/360;
                 [vfEstimate] = circ_ksdensity(vfObservations, vfPDFSamples, vfDomain, fSigma);
                 
                 currPdist=vfEstimate/sum(vfEstimate(:));
               
            else
                currPdist=circSmooth(pxySmooth(ti,:),3);
            end
            %currPdist=circSmooth(pxy(ti,:),5);
            %currPdist=circSmooth(pxy(ti,:),5);

            currPdist=currPdist/sum(currPdist);


              [~,maxID]=max(currPdist);
        %    currPdist=circshift(currPdist,-maxID);

           if(useKappaVarMeasure)
               currDistrWidth=sqrt(currTimeBinKappa);
           else
               upsampleFact=30;
               [phaseBinCentersPadded,currPdistPadded,padLength]=getCircPaddedDistr(phaseBinCenters,currPdist,upsampleFact);
         
         
           %[currDistrWidth,risingIntersectPhase,risingIntersectProb,descendingIntersectPhase,descendingIntersectProb]...
            %   =getFullWidthHalfMaxCircDistr(phaseBinCenters,currPdist);
           
           [currDistrWidth,risingIntersectPhase,risingIntersectProb,descendingIntersectPhase,descendingIntersectProb]...
               =getFullWidthFracMaxCircDistr(phaseBinCenters,currPdist,0.5);
           end
           
           
           if(currDistrWidth>360)
               continue
           end
        if(showPlots)
        subplot(1,3,1)
    plot(currFieldInFieldTimes,currFieldSpikePhases,'k.','MarkerSize',25)
       
    subplot(1,3,2)
         if(useKappaVarMeasure)
               %currDistrWidth=sqrt(currTimeBinKappa);
               plot(phaseBinCenters,currPdist,'-','Color',timeBinColors(ti,:),'LineWidth',3)
               xlim([0 360])
           else

            plot(phaseBinCentersPadded,currPdistPadded,'-','Color',timeBinColors(ti,:),'LineWidth',3)

            hold on

            plot(risingIntersectPhase,risingIntersectProb,'.','Color',timeBinColors(ti,:),'MarkerSize',30)
            plot(descendingIntersectPhase,descendingIntersectProb,'.','Color',timeBinColors(ti,:),'MarkerSize',30)

            plot([risingIntersectPhase descendingIntersectPhase], [risingIntersectProb descendingIntersectProb],'Color',timeBinColors(ti,:),'LineWidth',3)

         end


            %circMeanPhasePerSpeedBin(ti)=circMeanDeg(


            hold on
           subplot(1,3,3)
            
            if(takeSqrt)
                if(currDistrWidth>=baselineWidth)
                    currDistrWidthSqrt=sqrt(currDistrWidth-baselineWidth);
                else
                    currDistrWidthSqrt=0;
                end
      
                plot(timeBinCenters(ti),currDistrWidthSqrt,'.','Color',timeBinColors(ti,:),'MarkerSize',80)
        
            else
                
                plot(timeBinCenters(ti),currDistrWidth,'.','Color',timeBinColors(ti,:),'MarkerSize',80)
   
            end
            hold on
           
            %ylim([30 270])
            if(useKappaVarMeasure)
                ylim([0 3])
            else
                if(takeSqrt)
                    ylim([0 20])
                else
                %ylim([30 330])
                ylim([0 300])
                end
            end
            
            
    end
             if(takeSqrt)
                 distrWidth(ti)=currDistrWidthSqrt;
             else
                 distrWidth(ti)=currDistrWidth;
             end
            maxProbPhasePerTimeBin(ti)=phaseBinCenters(maxID);
    end
    
    notNaNIdxes=~isnan(distrWidth);
        timeBinCentersNonNan=timeBinCenters(notNaNIdxes);
        distrWidthNonNan=distrWidth(notNaNIdxes);
        
        currFieldAvgPhaseWidth=mean(distrWidthNonNan);
        %currFieldAvgPhaseWidth=0;
       %currFieldAvgPhaseWidth=nanmin(distrWidthNonNan);
        
        [m,b,R]=getLinearFit(timeBinCentersNonNan,distrWidthNonNan);
            y0=b;
            yf=m+b;
            
        %currFieldAvgPhaseWidth=b;
        
        predictedVals=interp1([0 1],[y0 yf],timeBinCentersNonNan);
        actualVals=distrWidthNonNan;

        signedLinErrorsRealMinusPredicted=(actualVals(:)-predictedVals(:));
        
       

        %timeBinCenters
        
        allFieldTimeBinCenters=[allFieldTimeBinCenters;timeBinCentersNonNan(:)];
        allFieldDistrWidths=[allFieldDistrWidths;distrWidthNonNan(:)];      
        allFieldLinearModelErrs=[allFieldLinearModelErrs;signedLinErrorsRealMinusPredicted(:)];
        
        allFieldStatIDs=[allFieldStatIDs;repelem(numFieldsUsed,length(distrWidthNonNan),1)];

    if(showPlots)
        
         subplot(1,3,1)
  
    xlabel('Time in field (frac)')
    ylabel('Theta phase (deg)')
    xlim([0 1])
    ylim([0 360])
        title(sprintf('Theta phase vs Time in traversal, field %d (n=%d spikes)',fi,length(currFieldSpikePhases)))

    subplot(1,3,2)
    colormap(gca,jet)
    cb=colorbar
    ylabel(cb,'Time in field (frac)')
    xlabel('Theta phase (deg)')

    ylabel('Probability')
    %xlim([0 360])
    xlim([-20 380])
    caxis([0 1])

    title(sprintf('Theta phase distribution vs Time in traversal'))


    %daspect([1 360 1])
    %caxis([0 prctile(pxySmooth(:),97.5)])
    %maxFig
    subplot(1,3,3)
    
    if(useKappaVarMeasure)
        title({sprintf('Distribution width vs time, R=%.2f',R),sprintf('slope=%.2f kappa/field',m)})
    else
       
       title({sprintf('Distribution width vs time, R=%.2f',R),sprintf('slope=%.2f deg/field',m)})
    end
    
        hold on
    x0=0;
    y0=b;
    xf=1;
    yf=m+b;
    plot([x0 xf], [y0 yf],'k--','LineWidth',5)
    
     xlabel('Time in field (frac)')
      if(useKappaVarMeasure)
          ylabel('sqrt(Kappa)')
      else
            ylabel('Theta phase distribution width (deg)')
      end
      

     setFigFontTo(16)
    maxFig
    saveas(gcf,fullfile(saveFieldLogCompressPhaseVarDir,sprintf('timeVsThetaPhaseJointDistr_Field%d.png',fi)))
    close all

    
    end
    
    allFieldTimeWidthSlopes=[allFieldTimeWidthSlopes; m];
    allFieldTimeWidthRs=[allFieldTimeWidthRs; R];
    allFieldTimeWidthAvgs=[allFieldTimeWidthAvgs; currFieldAvgPhaseWidth];

end

%%
numFieldsUsed

save(saveDataPath,'allFieldStatIDs', 'allFieldTimeWidthAvgs', 'allFieldTimeWidthSlopes','allFieldTimeWidthRs','allFieldTimeBinCenters','allFieldDistrWidths','minNumSpikesPerTimeBin','timeBinCenters','allFieldLinearModelErrs')






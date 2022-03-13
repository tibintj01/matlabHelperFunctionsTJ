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
numTimeBins=20;
%numTimeBins=10;

timeEdges=linspace(0,1,numTimeBins+1);
%timeEdges=linspace(0,1.5,numTimeBins+1);
timeBinCenters=edgesToBins(timeEdges);
timeBinWidth=median(diff(timeBinCenters));


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
numPhaseBins=20;
%phaseEdges=linspace(0,360,numPhaseBins+1);

phaseEdges=linspace(40,360,numPhaseBins+1);

phaseBinCenters=edgesToBins(phaseEdges);

phaseBinWidth=median(diff(phaseBinCenters))

phaseBinColors=jet(numPhaseBins)
maxProbTimePerPhaseBin=NaN(numTimeBins,1);
distrWidth=NaN(numPhaseBins,1);

minNumSpikesForField=600;
minNumSpikesForField=700;
minNumSpikesForField=300;


minNumSpikesPerPhaseBin=10;
minNumSpikesPerTimeBin=30;
%minNumSpikesForField=3000;

numFieldsUsed=0;

showPlots=0;
showPlots=1;

allFieldTimeWidthSlopes=[];
allFieldTimeWidthRs=[];
allFieldTimeBinCenters=[];
allFieldTimeDistrWidths=[];
allFieldLinearModelErrs=[];

for fi=1:totalFieldCount
    
    currFieldSpikePhases=spikePhasesPerField{fi};
    currFieldInFieldTimes=spikeTimeFracInFieldPerField{fi};
    
    if(length(currFieldSpikePhases)<minNumSpikesForField)
        continue
    else
          
        numFieldsUsed=numFieldsUsed+1;
        %continue
    end
    
    

    [pxySmooth,pxy]=getJointDistr(currFieldInFieldTimes,currFieldSpikePhases,timeEdges,phaseEdges);

    if(showPlots)
    figure;
    title(sprintf('Phase distribution vs time in field, field %d',fi))
    xlabel('Time in field (frac)')
    ylabel('Theta phase (deg)')
    end
    
    

    for phi=1:numPhaseBins
        m=NaN;
        R=NaN;
        b=NaN;
        
        %get kappa for current bin
         currPhaseBinStart=phaseBinCenters(phi)-phaseBinWidth/2;
         currPhaseBinEnd=phaseBinCenters(phi)+phaseBinWidth/2;
            
            inBinIdxes=currFieldSpikePhases>=currPhaseBinStart & currFieldSpikePhases<currPhaseBinEnd;
            
          
           inBinTimes=currFieldInFieldTimes(inBinIdxes);
           notNaNIdxesTimes=~isnan(inBinTimes) & inBinTimes<=1 & inBinTimes>=0;
           
           if(sum(notNaNIdxesTimes)<=minNumSpikesPerPhaseBin )
               continue
           end
           
    
           inBinTimes=inBinTimes(notNaNIdxesTimes);
            
            
            %currTimeBinKappa=circ_kappa(ang2rad(inBinTimes(notNaNIdxesTimes)));
            
            if(useKernelEstDistr)
                kernelDistrTimeObj = fitdist(inBinTimes(:),'Kernel','BandWidth',0.1); %gaussian kernel estimation
                %phaseAxis=linspace(0,360,361);
                %phaseAxis=linspace(0,360,numPhaseBins+1);
                currPdist=pdf(kernelDistrTimeObj,timeBinCenters);
            else
                currPdist=circSmooth(pxySmooth(phi,:),3);
            end
            %currPdist=circSmooth(pxy(ti,:),5);
            %currPdist=circSmooth(pxy(ti,:),5);

            currPdist=currPdist/sum(currPdist);


              [~,maxID]=max(currPdist);
        %    currPdist=circshift(currPdist,-maxID);

        
        
         
           
           if(useKappaVarMeasure)
               currDistrWidth=sqrt(currTimeBinKappa);
           else
               %[phaseBinCentersPadded,currPdistPadded,padLength]=getCircPaddedDistr(phaseBinCenters,currPdist,10);
         
         
               [currDistrWidth,risingIntersectTime,risingIntersectProb,descendingIntersectTime,descendingIntersectProb]...
                   =getFullWidthHalfMaxCircDistr(timeBinCenters*360,currPdist);

               currDistrWidth=currDistrWidth/360;

                risingIntersectTime=risingIntersectTime/360;
               descendingIntersectTime=descendingIntersectTime/360;
           end
           
           
           if(currDistrWidth>360)
               continue
           end
           
        if(showPlots)
        subplot(1,2,1)

       
         if(useKappaVarMeasure)
               %currDistrWidth=sqrt(currTimeBinKappa);
               plot(phaseBinCenters,currPdist,'-','Color',timeBinColors(phi,:),'LineWidth',3)
               xlim([0 360])
           else

             plot(timeBinCenters,currPdist,'-','Color',phaseBinColors(phi,:),'LineWidth',3)

            hold on

            plot(risingIntersectTime,risingIntersectProb,'.','Color',phaseBinColors(phi,:),'MarkerSize',30)
            plot(descendingIntersectTime,descendingIntersectProb,'.','Color',phaseBinColors(phi,:),'MarkerSize',30)

            plot([risingIntersectTime descendingIntersectTime], [risingIntersectProb descendingIntersectProb],'Color',phaseBinColors(phi,:),'LineWidth',3)


         end


            %circMeanPhasePerSpeedBin(ti)=circMeanDeg(


            hold on
            subplot(1,2,2)
            
           
                
           plot(phaseBinCenters(phi),currDistrWidth,'.','Color',phaseBinColors(phi,:),'MarkerSize',80)
   
   
            hold on
           
            %ylim()
            
           
        end
          
             distrWidth(phi)=currDistrWidth;
            
            maxProbTimePerPhaseBin(phi)=timeBinCenters(maxID);
    end
    
    notNaNIdxes=~isnan(distrWidth) & distrWidth<=1 & distrWidth>0;
        phaseBinCentersNonNan=phaseBinCenters(notNaNIdxes);
        distrWidthNonNan=distrWidth(notNaNIdxes);
        [m,b,R]=getLinearFit(phaseBinCentersNonNan,distrWidthNonNan);
            y0=b;
            yf=m*360+b;
        
        predictedVals=interp1([0 360],[y0 yf],phaseBinCentersNonNan);
        actualVals=distrWidthNonNan;

        signedLinErrorsRealMinusPredicted=(actualVals(:)-predictedVals(:));
        
        
     

        %timeBinCenters
        
        allFieldTimeBinCenters=[allFieldTimeBinCenters;phaseBinCentersNonNan(:)];
        allFieldTimeDistrWidths=[allFieldTimeDistrWidths;distrWidthNonNan(:)];
        
        allFieldLinearModelErrs=[allFieldLinearModelErrs;signedLinErrorsRealMinusPredicted(:)];

    if(showPlots)
    subplot(1,2,1)
    colormap(gca,jet)
    cb=colorbar
    ylabel(cb,'Theta phase (deg)')
    xlabel('Time in field (frac)')

    ylabel('Probability')
    %xlim([0 360])
    xlim([0 1])
    caxis([0 360])

    title(sprintf('Time distribution vs theta phase, field %d (n=%d spikes)',fi,length(currFieldSpikePhases)))


    %daspect([1 360 1])
    %caxis([0 prctile(pxySmooth(:),97.5)])
    %maxFig
    subplot(1,2,2)
  
    title({sprintf('Time distribution width vs phase, R=%.2f',R),sprintf('slope=%.2f time/deg',m)})

   
    
        hold on
    x0=0;
    y0=b;
    xf=360;
    yf=m*360+b;
    plot([x0 xf], [y0 yf],'k--','LineWidth',5)
    
     xlabel('Theta phase (deg)')
     
    ylabel('Time distribution width (frac)')
    ylim([0 1])
    xlim([0 360])
 
      
     setFigFontTo(18)
    maxFig
    saveas(gcf,fullfile(saveFieldLogCompressPhaseVarDir,sprintf('timeVsThetaPhaseJointDistr_Field%d.png',fi)))
    close all

    
    end
    
    allFieldTimeWidthSlopes=[allFieldTimeWidthSlopes; m];
    allFieldTimeWidthRs=[allFieldTimeWidthRs; R];
  

end

%%
numFieldsUsed

save(saveDataPath,'allFieldTimeWidthSlopes','allFieldTimeWidthRs','allFieldTimeBinCenters','allFieldTimeDistrWidths','minNumSpikesPerTimeBin','phaseBinCenters','allFieldLinearModelErrs')






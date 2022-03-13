close all;clear all; clc

tic
disp('loading per field data...')
dataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
toc
totalFieldCount=dataPerField.totalFieldCount;
speedVsPhaseInfoPerField=dataPerField.speedVsPhaseInfoPerField;

typicalSpeedBinEdges=[0.1 0.2 0.3 0.4 0.5 0.6];

typicalSpeedBinEdges=[0.1 0.3 0.4 0.5 0.6];

typicalSpeedBinEdges=[0.3 0.5];

%typicalSpeedBinEdges=[0 0.2 0.4 0.6 0.8];

typicalSpeedBinEdges=[0 1];

%typicalSpeedBinEdges=[0.2 0.4 0.5 0.6 0.7 0.8];

typicalSpeedBinCenters=edgesToBins(typicalSpeedBinEdges);

numCharSpeedBins=length(typicalSpeedBinCenters);

showPlots=1;

maxDispTime=1;

fH=figure
allLowSpeedSpikeTimeFracs=cell(numCharSpeedBins,1);
allLowSpeedSpikePhases=cell(numCharSpeedBins,1);
allHighSpeedSpikeTimeFracs=cell(numCharSpeedBins,1);
allHighSpeedSpikePhases=cell(numCharSpeedBins,1);

timeInFieldEdges=linspace(0,1.2,31);
phaseEdges=linspace(0,360,101);

setTightSubplots_Medium
       
for fi=1:totalFieldCount
    fi
    currFieldSpeedVsPhaseInfo=speedVsPhaseInfoPerField{fi};
    
    if(isempty(currFieldSpeedVsPhaseInfo))
        continue
    end
    
    charSpeed=currFieldSpeedVsPhaseInfo.characteristicSpeedMperSec;
    meanOffset=currFieldSpeedVsPhaseInfo.allEarlyFieldMeanOffset;
    
    highMinusLowPhaseDiffPerTimeFrac=currFieldSpeedVsPhaseInfo.highMinusLowPhaseDiffPerTimeFrac;
    uniqueTimeFracs=currFieldSpeedVsPhaseInfo.uniqueTimeFracs;
    
    currFieldLowSpeedSpikeTimeFracs=currFieldSpeedVsPhaseInfo.currFieldLowSpeedSpikeTimeFracs;
    currFieldLowSpeedSpikePhases=currFieldSpeedVsPhaseInfo.currFieldLowSpeedSpikePhases;
    currFieldHighSpeedSpikeTimeFracs=currFieldSpeedVsPhaseInfo.currFieldHighSpeedSpikeTimeFracs;
    currFieldHighSpeedSpikePhases=currFieldSpeedVsPhaseInfo.currFieldHighSpeedSpikePhases;
    
    charSpeedBin=discretize(charSpeed,typicalSpeedBinEdges);
    
    if(isnan(charSpeedBin))
        %continue
    end
    
    allLowSpeedSpikeTimeFracs{charSpeedBin}=[allLowSpeedSpikeTimeFracs{charSpeedBin};currFieldLowSpeedSpikeTimeFracs(:)];
     allLowSpeedSpikePhases{charSpeedBin}=[allLowSpeedSpikePhases{charSpeedBin};currFieldLowSpeedSpikePhases(:)];
      allHighSpeedSpikeTimeFracs{charSpeedBin}=[allHighSpeedSpikeTimeFracs{charSpeedBin};currFieldHighSpeedSpikeTimeFracs(:)];
       allHighSpeedSpikePhases{charSpeedBin}=[allHighSpeedSpikePhases{charSpeedBin};currFieldHighSpeedSpikePhases(:)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(showPlots)% && (mod(fi,100)==0 || fi==totalFieldCount))
        
        %{
        subplot(numCharSpeedBins,2,2*(charSpeedBin-1)+1)
        %plot(currFieldLowSpeedSpikeTimeFracs,currFieldLowSpeedSpikePhases,'b.')
        %getJointDistr(x,y,xEdges,yEdges,fH,withSmooth,cmapName)
        getJointDistr(allLowSpeedSpikeTimeFracs{charSpeedBin},allLowSpeedSpikePhases{charSpeedBin},timeInFieldEdges,phaseEdges,fH,0,jet)
        %getJointDistrGivenX(allLowSpeedSpikeTimeFracs,allLowSpeedSpikePhases,timeInFieldEdges,phaseEdges,fH,0,jet)

        title(sprintf('%.2f < Char speed < %.2f, low speeds',typicalSpeedBinEdges(charSpeedBin),typicalSpeedBinEdges(charSpeedBin+1)))
        hold on
        xlim([0 maxDispTime])
        ylim([0 360])
          axis square
        
        subplot(numCharSpeedBins,2,2*(charSpeedBin-1)+2)
        %plot(currFieldHighSpeedSpikeTimeFracs,currFieldHighSpeedSpikePhases,'r.')
        %getJointDistrGivenX(allHighSpeedSpikeTimeFracs,allHighSpeedSpikePhases,timeInFieldEdges,phaseEdges,fH,0,jet)
        getJointDistr(allHighSpeedSpikeTimeFracs{charSpeedBin},allHighSpeedSpikePhases{charSpeedBin},timeInFieldEdges,phaseEdges,fH,0,jet)

        
        title(sprintf('%.2f < Char speed < %.2f field, high speeds',typicalSpeedBinEdges(charSpeedBin),typicalSpeedBinEdges(charSpeedBin+1)))
        hold on
        xlim([0 maxDispTime])
        ylim([0 360])
          axis square
        
       
        %}
        
        
        
      
        
        plot(uniqueTimeFracs,highMinusLowPhaseDiffPerTimeFrac,'k.')
        xlabel('Time in field (frac)')
        ylabel('Avg high speed phase - avg low speed phase (deg)')
        title('High speed phase difference across field')
        xlim([0 1.2])
        ylim([-180 180])
        
        %{
        plot(charSpeed,meanOffset,'k.')
        
        xlabel('Characteristic speed of field (m/s)')
        ylabel('circ mean early field phase (deg)')
        title('Characteristic speed vs offset confound when pooling fields')

        %}
         hold on
        drawnow
        
        
    end
    
    
    
end

[timeFracBinCenters,lowSpeedPhaseAvgValues]= ...
            getBinnedCircAverages(allLowSpeedSpikeTimeFracs{charSpeedBin},allLowSpeedSpikePhases{charSpeedBin},timeInFieldEdges)
        
         [timeFracBinCenters,highSpeedPhaseAvgValues]= ...
            getBinnedCircAverages(allHighSpeedSpikeTimeFracs{charSpeedBin},allHighSpeedSpikePhases{charSpeedBin},timeInFieldEdges)
        
        figure;
        
        plot(timeFracBinCenters,lowSpeedPhaseAvgValues,'b')
        hold on
         plot(timeFracBinCenters,lowSpeedPhaseAvgValues,'b.','MarkerSize',20)

         plot(timeFracBinCenters,highSpeedPhaseAvgValues,'r')
        hold on
         plot(timeFracBinCenters,highSpeedPhaseAvgValues,'r.','MarkerSize',20)
         
         xlabel('Time in field (frac)')
         ylabel('Avg theta phase (deg)')
         ylim([0 360])
         xlim([0 1])
    


setFigFontTo(14)
%saveas(gcf,'charSpeedControlledHighVsLowSpeedPhasePrecession.png')
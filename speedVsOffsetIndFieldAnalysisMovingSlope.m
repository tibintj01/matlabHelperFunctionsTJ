close all;clear all; clc

tic
disp('loading per field data...')
dataPerField=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');
toc
totalFieldCount=dataPerField.totalFieldCount;
speedVsPhaseInfoPerField=dataPerField.speedVsPhaseInfoPerField;

typicalSpeedBinEdges=[0.1 0.2 0.3 0.4 0.5 0.6];


typicalSpeedBinCenters=edgesToBins(typicalSpeedBinEdges);

numCharSpeedBins=length(typicalSpeedBinCenters);

showPlots=1;

maxDispTime=1;
%maxDispTime=1.2;

fH=figure
%{
allLowSpeedSpikeTimeFracs=cell(numCharSpeedBins,1);
allLowSpeedSpikePhases=cell(numCharSpeedBins,1);
allHighSpeedSpikeTimeFracs=cell(numCharSpeedBins,1);
allHighSpeedSpikePhases=cell(numCharSpeedBins,1);
%}

allLowSpeedSpikeTimeFracs=[];
allLowSpeedSpikePhases=[];
allHighSpeedSpikeTimeFracs=[];
allHighSpeedSpikePhases=[];

timeInFieldEdges=linspace(0,1.2,31);
phaseEdges=linspace(0,360,101);

setTightSubplots_Medium
       
%allTimeFracVsPhaseDiffRperField=NaN(totalFieldCount,1);
    allTimeFracVsPhaseDiffRperField=dataPerField.allTimeFracVsPhaseDiffRperField;
timeFieldCount=0;
for fi=1:totalFieldCount
    fi
    currFieldSpeedVsPhaseInfo=speedVsPhaseInfoPerField{fi};
    
    if(isempty(currFieldSpeedVsPhaseInfo))
        continue
    end
    
    charSpeed=currFieldSpeedVsPhaseInfo.characteristicSpeedMperSec;
    meanOffset=currFieldSpeedVsPhaseInfo.allEarlyFieldMeanOffset;
    
    allTimeFracVsPhaseDiffR=allTimeFracVsPhaseDiffRperField(fi);
    
    highMinusLowPhaseDiffPerTimeFrac=currFieldSpeedVsPhaseInfo.highMinusLowPhaseDiffPerTimeFrac;
    uniqueTimeFracs=currFieldSpeedVsPhaseInfo.uniqueTimeFracs;
    
    %if(allTimeFracVsPhaseDiffR<0 || nanmean(highMinusLowPhaseDiffPerTimeFrac) >0)
    if(allTimeFracVsPhaseDiffR<0 || nanmean(highMinusLowPhaseDiffPerTimeFrac) >0)
    %if(allTimeFracVsPhaseDiffR<-0.15 || allTimeFracVsPhaseDiffR>0.15 )
        continue
    end
    
    timeFieldCount=timeFieldCount+1;
    
    
    currFieldLowSpeedSpikeTimeFracs=currFieldSpeedVsPhaseInfo.currFieldLowSpeedSpikeTimeFracs;
    currFieldLowSpeedSpikePhases=currFieldSpeedVsPhaseInfo.currFieldLowSpeedSpikePhases;
    currFieldHighSpeedSpikeTimeFracs=currFieldSpeedVsPhaseInfo.currFieldHighSpeedSpikeTimeFracs;
    currFieldHighSpeedSpikePhases=currFieldSpeedVsPhaseInfo.currFieldHighSpeedSpikePhases;
    
    allLowSpeedSpikeTimeFracs=[allLowSpeedSpikeTimeFracs; currFieldLowSpeedSpikeTimeFracs(:)];
    allLowSpeedSpikePhases=[allLowSpeedSpikePhases; currFieldLowSpeedSpikePhases(:)];
    allHighSpeedSpikeTimeFracs=[allHighSpeedSpikeTimeFracs; currFieldHighSpeedSpikeTimeFracs(:)];
    allHighSpeedSpikePhases=[allHighSpeedSpikePhases; currFieldHighSpeedSpikePhases(:)];
    


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(showPlots)
      
        
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
        
        
        
        
      %[m,b,R,p]=getLinearFit(uniqueTimeFracs,highMinusLowPhaseDiffPerTimeFrac);
      
      
        %allTimeFracVsPhaseDiffRperField(fi)=R;
      %{
        plot(uniqueTimeFracs,highMinusLowPhaseDiffPerTimeFrac,'k.')
        xlabel('Time in field (frac)')
        ylabel('Avg high speed phase - avg low speed phase (deg)')
        title('High speed phase difference across field')
        xlim([0 1.2])
        ylim([-180 180])
      %}
        
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

  subplot(1,2,1)
         getJointDistr(allLowSpeedSpikeTimeFracs,allLowSpeedSpikePhases,timeInFieldEdges,phaseEdges,fH,0,jet)
         title(sprintf('Time phase-coding fields, low speeds (n=%d fields)',timeFieldCount))
         xlim([0 maxDispTime])
        ylim([0 360])
          axis square
           caxis([0 1.5e-3])
        subplot(1,2,2)
         getJointDistr(allHighSpeedSpikeTimeFracs,allHighSpeedSpikePhases,timeInFieldEdges,phaseEdges,fH,0,jet)
       
        title(sprintf('Time phase-coding fields, high speeds (n=%d fields)',timeFieldCount))
        xlim([0 maxDispTime])
        ylim([0 360])
          axis square
          caxis([0 1.6e-3])
setFigFontTo(14)
saveas(gcf,'highVsLowSpeedPhasePrecessionTimeCells_JointDistr.png')
%save('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat','allTimeFracVsPhaseDiffRperField','-append')
%%
timeInFieldEdges=linspace(0,maxDispTime,21);
[timeFracBinCenters,lowSpeedPhaseAvgValues]= ...
            getBinnedCircAverages(allLowSpeedSpikeTimeFracs,allLowSpeedSpikePhases,timeInFieldEdges)
        
         [timeFracBinCenters,highSpeedPhaseAvgValues]= ...
            getBinnedCircAverages(allHighSpeedSpikeTimeFracs,allHighSpeedSpikePhases,timeInFieldEdges)
        
        figure;
        
        pB=plot(timeFracBinCenters,lowSpeedPhaseAvgValues,'b')
        hold on
         plot(timeFracBinCenters,lowSpeedPhaseAvgValues,'b.','MarkerSize',20)

         pR=plot(timeFracBinCenters,highSpeedPhaseAvgValues,'r')
        hold on
         plot(timeFracBinCenters,highSpeedPhaseAvgValues,'r.','MarkerSize',20)
         
         xlabel('Time in field (frac)')
         ylabel('Avg theta phase (deg)')
         ylim([0 360])
         xlim([0 maxDispTime])
    
         legend([pB pR],{'low speeds', 'high speeds'})
         
          title(sprintf('Time phase-coding only (n=%d fields)',timeFieldCount))


setFigFontTo(14)
saveas(gcf,'highVsLowSpeedPhasePrecessionTimeCells_Avg.png')
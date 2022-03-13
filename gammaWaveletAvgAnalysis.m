clear all;close all; clc
allFieldData=load('allUnitInFieldGoodSpeedSpikePhasesAndTimes.mat');

totalFieldCount=allFieldData.totalFieldCount;
wtCommonGammaFreqAxis=allFieldData.wtCommonGammaFreqAxis;
wtCommonThetaPhaseAxis=allFieldData.wtCommonThetaPhaseAxis;

periSpikeGammaWTacrossThetaPerField=allFieldData.periSpikeGammaWTacrossThetaPerField_JustAroundSpikes;
fH=figure;
fHist=figure;

allFieldGammaWTzScoreAvg=zeros(size(periSpikeGammaWTacrossThetaPerField{1}));

for fi=1:totalFieldCount
    currFieldGammaWT=periSpikeGammaWTacrossThetaPerField{fi};
    
    meanPower=nanmean(currFieldGammaWT(:));
    stdPower=nanstd(currFieldGammaWT(:));
    
    currFieldGammaWTzScore=(currFieldGammaWT-meanPower)/stdPower;
    
    meanTimeCourse=nanmean(currFieldGammaWTzScore,1);
    
    
    currFieldGammaWTzScore=currFieldGammaWTzScore-meanTimeCourse;
    
    meanSpectrumAcrossPhase=nanmean(currFieldGammaWTzScore,2);
    currFieldGammaWTzScore=currFieldGammaWTzScore-meanSpectrumAcrossPhase;
    
    currFieldGammaWTzScore(isnan(currFieldGammaWTzScore))=0;
    
    
    
    allFieldGammaWTzScoreAvg=allFieldGammaWTzScoreAvg+currFieldGammaWTzScore;
    
    
    %allFieldGammaWTzScoreAvg=allFieldGammaWTzScoreAvg/fi;
    
    figure(fH)
    %subplot(2,1,1)
    omarPcolor(wtCommonThetaPhaseAxis,wtCommonGammaFreqAxis,allFieldGammaWTzScoreAvg/fi,fH)
        %omarPcolor(wtCommonThetaPhaseAxis,wtCommonGammaFreqAxis,currFieldGammaWTzScore,fH)

    colormap(jet)
    colorbar
    
    %caxis([-3 3])
    
    %pause(0.4)
    
    %{
    subplot(2,1,2)
    omarPcolor(wtCommonThetaPhaseAxis,wtCommonGammaFreqAxis,allFieldGammaWTzScoreAvg,fH)
    colormap(jet)
    colorbar
    drawnow
    %}
    %{
    figure(fH)
    subplot(10,10,fi)
    omarPcolor(wtCommonThetaPhaseAxis,wtCommonGammaFreqAxis,currFieldGammaWTzScore,fH)
    colormap(jet)
    colorbar
    
    figure(fHist)
    subplot(10,10,fi)
    histogram(currFieldGammaWTzScore)
    
    if(fi==100)
        
      disp('')
      
    end
    %}
end

title(sprintf('Avg peri-spike wavelet transform within theta cycles, n=%d fields',totalFieldCount))
xlabel('Theta phase (deg)')

ylabel('Frequency (Hz)')
cb=colorbar;
ylabel(cb,'Power (z-score relative to non-spike times)')
setFigFontTo(18)
%%
%binned average plot
phaseColors=jet(length(wtCommonThetaPhaseAxis));
figure
comGammaFreqPerThetaPhase=NaN(size(wtCommonThetaPhaseAxis));

minGammaFreq=40;
%maxGammaFreq=75;
maxGammaFreq=120;

for ti=1:1:length(wtCommonThetaPhaseAxis)
    
    currPhaseGammaSpec=allFieldGammaWTzScoreAvg(:,ti)/fi;
    currPhaseGammaSpec(wtCommonGammaFreqAxis<minGammaFreq | wtCommonGammaFreqAxis>maxGammaFreq)=NaN;
    %plot(wtCommonGammaFreqAxis,currPhaseGammaSpec,'Color',phaseColors(ti,:))
    %hold on
    
    currPhaseCOMGammafreqIdx=getCOM(currPhaseGammaSpec-min(currPhaseGammaSpec));
    
    currPhaseCOMGammafreq=interp1(1:length(wtCommonGammaFreqAxis), wtCommonGammaFreqAxis,currPhaseCOMGammafreqIdx);
    comGammaFreqPerThetaPhase(ti)=currPhaseCOMGammafreq;
    
end
saveas(gcf,'avgPerSpikeWaveletTransformWithinThetaCycles.png')

figure; 
plot([wtCommonThetaPhaseAxis-360 wtCommonThetaPhaseAxis wtCommonThetaPhaseAxis+360],[comGammaFreqPerThetaPhase comGammaFreqPerThetaPhase comGammaFreqPerThetaPhase],'k-','LineWidth',5)
hold on
vline(360,'k--')
vline(0,'k--')

xlabel('Theta phase (deg)')
ylabel('Center of mass slow-mid gamma frequency (Hz)')

title(sprintf('Slow to mid gamma frequency (%d-%d Hz) vs theta phase',minGammaFreq,maxGammaFreq))
setFigFontTo(18)

buffer=360;
xlim([-buffer 360+buffer])

saveas(gcf,'comSlowMidGammaFreqVsThetaPhase.png')

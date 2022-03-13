
singleCycleType='alpha6-20Hz';
singleCycleType='alpha6-40Hz';
%singleCycleType='alpha8_2-13Hz';

singleCycleMatDirPath=sprintf('/nfs/turbo/lsa-ojahmed/processedHumanData/%s/sessionID-%d/singleCycleProperties-MatFiles/%s',ptDir,sessionNum,singleCycleType);
saveDir=fullfile(singleCycleMatDirPath,[singleCycleType '-cycleMatrix'])

cycleFeatureMathPath=sprintf('/nfs/turbo/lsa-ojahmed/processedHumanData/MG49/sessionID-3/singleCycleProperties-MatFiles/%s/AlphaSingleCycleBroadbandStats-Ch06.mat',singleCycleType);

data=load(cycleFeatureMathPath);

durations=data.asymData.fastStartEndDuration;

figure
subplot(2,2,1)
histogram(durations)
title({singleCycleType,'MG49-3, Ch 06 Durations'})
xlabel('Duration (sec)')
xlim([0 0.4])

troughPhases=data.asymData.fastStartMidDuration./durations*360;
subplot(2,2,2)
histogram(troughPhases)
title({singleCycleType,'Trough Phases'})
xlabel('Phase of midpoint (degrees)')
xlim([0 360])

amps=data.asymData.fastStartMidAmp;
subplot(2,2,3)
histogram(amps)
title({singleCycleType,'Start to Mid Cycle amplitudes'})
xlabel('Amplitude')
xlim([0 30])

slopeRatios=data.asymData.fastSlopeRatio;
subplot(2,2,4)
histogram(slopeRatios)
title({singleCycleType,'Slope ratio'})
xlabel('Finish to start slopes ratio')
xlim([0 10])

saveas(gcf,fullfile(saveDir,sprintf('cycleFeatureDistributions%s.tif',singleCycleType)))


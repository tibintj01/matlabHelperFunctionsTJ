function [spikeWaveforms] = getWaveformsFromLFP(lfpTimeSeries,spikeTimes,deadSpace)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%
%
%Input:
%
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-09-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lfp=lfpTimeSeries.lfp;
Fs=lfpTimeSeries.Fs;

%samplesBeforeCrossing=9;
%samplesAfterCrossing=38;
samplesBeforeCrossing=21;
samplesAfterCrossing=38;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing getRawWaveformsFromLFP fxn output.........')


        %deadSpace=500; %16.67ms on either side (5 cycles of 300 Hz to avoid edge effects) - enough?
        realFirstIdx=deadSpace+1;
        realLastIdx=realFirstIdx+samplesBeforeCrossing+1+samplesAfterCrossing-1;

        spikeWaveforms=zeros(length(spikeTimes),deadSpace+48+deadSpace);
        %spikeWaveforms=NaN(length(spikeTimes),deadSpace+48+deadSpace);

        tic
        disp('getting filt segment for each spike....')
        numWithinTimeSpikes=0;
        for spikeNum=1:length(spikeTimes)
                spikeTime=spikeTimes(spikeNum);
                startTime=spikeTime-samplesBeforeCrossing/Fs;
                %startTimeMinus=max(0,startTime-2);

                endTime=spikeTime+samplesAfterCrossing/Fs;
                %endTimePlus=min(endTime+2,length(filtLFP)/Fs);
                %filtSpikeWaveform=getRawSegment_MG49_sess3(startTime,endTime,chanNum);
                %filtSpikeWaveform=getRawSegmentFromConcat(startTime,endTime,chanNum);

                startSample=max(1,round(startTime*Fs)-deadSpace);
                endSample=min(length(lfp),round(endTime*Fs)+deadSpace);

                %startSampleMinus=round(startTimeMinus*Fs);
                %endSamplePlus=round(endTimePlus*F

                if(startSample > length(lfp) || endSample > length(lfp))
                        continue
                end
                %fds
                numWithinTimeSpikes=numWithinTimeSpikes+1;

                rawSpikeWaveform=lfp(startSample:endSample);
                rawSpikeZeroed=rawSpikeWaveform-median(rawSpikeWaveform(realFirstIdx:realLastIdx));
               %rawSpikeZeroed=rawSpikeWaveform;

                 spikeWaveforms(spikeNum,1:length(rawSpikeWaveform))=rawSpikeZeroed(:)';

                %filtSpikeWaveformsPlusMinus2msMat(spikeNum, 1:length(filtSpikeWaveformsPlusMinus2ms))=filtSpikeWaveformsPlusMinus2ms(:)';

        end
        toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikeWaveforms=spikeWaveforms';

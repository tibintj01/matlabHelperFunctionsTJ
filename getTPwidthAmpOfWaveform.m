function [spikeTPWidth,spikeTPAmp] = getTPwidthAmpOfWaveform(waveform)
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
%Created on 2018-10-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('computing getTPwidthAmpOfWaveform fxn output.........')

spikeLow=waveform;
numBins=length(waveform);

splineDt = .1; % XXX in bins (need to adjust based on sampling rate?)
 spikeHigh = spline(1:numBins, spikeLow, 1:splineDt:numBins);
    [spikeMinAmp spikeMinIndex] = min(spikeHigh); % should I restrict to first N bins? this should be a NEG value


spikeHighSmooth = fastsmooth(spikeHigh, 20, 3, 1); % NEW OMAR 2014 addition: Replace with causalSmooth?
    diffSpikeHigh = [0 diff(spikeHighSmooth)];

[ignoreMe spikeMinIndexSmooth] = min(spikeHighSmooth);
    firstDecreasingIndex = find(diffSpikeHigh(spikeMinIndexSmooth+1:end) <= 0, 1);

if isempty(firstDecreasingIndex) == 0
      spikeMaxIndex = firstDecreasingIndex + spikeMinIndexSmooth - 1; % XXX errorcheck this

	if (spikeMaxIndex - spikeMinIndex) > 10 % at least 1 original bin long
            spikeMaxAmp = spikeHigh(spikeMaxIndex);
            spikeTPAmp = spikeMaxAmp - spikeMinAmp;
            spikeTPWidth = (spikeMaxIndex-spikeMinIndex+1) * splineDt;
	else
	    spikeTPAmp = NaN;
	    spikeTPWidth = NaN;
	end
else
      spikeTPAmp = NaN;
      spikeTPWidth = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


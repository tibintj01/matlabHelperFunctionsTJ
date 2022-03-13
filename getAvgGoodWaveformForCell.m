function [avgGoodWaveform] = getGoodWaveforms(filePath)

	cellProp=load(filePath);

	%szStart=getSeizureStartTimes(filePath);
        %szStartTime=szStart(1);
        %goodSpikesB4sz=intersect(cellProp.goodSpikes,(find(cellProp.spikeTimes<szStartTime)));

        %goodWaveforms=cellProp.waveforms(:,goodSpikesB4sz);
	%goodWaveformsInterp=zeros(471,size(goodWaveforms,2));

        %numBins=size(goodWaveforms,1);
        %splineDt=0.1;
        %for i=1:size(goodWaveforms,2)
        %        wavInterped = spline(1:numBins, goodWaveforms(:,i), 1:splineDt:numBins);
        %        wavInterpedSmooth = fastsmooth(wavInterped, 20, 3, 1);
                %wavInterpedSmooth = wavInterped;
        %        goodWaveformsInterp(:,i)=wavInterpedSmooth;
        %end
	
	[goodWaveformsInterpTalign] = getGoodWaveforms(cellProp,filePath);

	avgGoodWaveform=nanmean(goodWaveformsInterpTalign,2);
	
	%if(size(avgGoodWaveform,1)<2)
	%	fds
	%end
end

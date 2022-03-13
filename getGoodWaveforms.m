function [goodWaveformsInterpTalign] = getGoodWaveforms(fileName)
	cellProp=load(fileName);
	szStart=getSeizureStartTimes(fileName);
        szStartTime=szStart(1);
        goodSpikesB4sz=intersect(cellProp.goodSpikes,find(cellProp.spikeTimes<szStartTime));

	minTimeBefore=0.010;
	minTimeAfter=0.005;
	isolatedSpikeIdxes=getIsolatedSpikeIdxes(cellProp.spikeTimes,minTimeBefore,minTimeAfter);	
	
	[allSpikeShifts,alignableSpikes]=getSpikeShifts(cellProp.spikeMinIndex);

	goodIsolatedSpikes=intersect(goodSpikesB4sz,isolatedSpikeIdxes);

	%goodTPwidths=cellProp.spikeTPWidth(goodIsolatedSpikes);
	%tpWidthFloor=prctile(goodTPwidths,90);
	%goodIsolatedSpikes=goodIsolatedSpikes(goodTPwidths>tpWidthFloor);
	goodMinAmps=cellProp.spikeMinAmp(goodIsolatedSpikes);
	ampPrc=30;
	ampCeiling=prctile(goodMinAmps,ampPrc);

	goodMaxAmps=cellProp.spikeMaxAmp(goodIsolatedSpikes);
	ampFloor=prctile(goodMaxAmps,100-ampPrc);
	
	%fds
	goodIsolatedSpikes=goodIsolatedSpikes(goodMinAmps<ampCeiling & goodMaxAmps>ampFloor);
	%tpWidths=cellProp.

	%goodIsolatedAlignableSpikes=intersect(goodIsolatedSpikes,alignableSpikes);
	%goodSpikeShifts=allSpikeShifts(goodIsolatedAlignableSpikes);
        %goodWaveforms=cellProp.waveforms(:,goodSpikesB4sz);
        goodWaveforms=cellProp.waveforms(:,goodIsolatedSpikes);
        %goodWaveforms=cellProp.waveforms(:,goodIsolatedAlignableSpikes);


        numBins=size(goodWaveforms,1);
        splineDt=0.1;
	
	goodWaveformsInterp=NaN(length(1:splineDt:numBins),size(goodWaveforms,2));
        for i=1:size(goodWaveforms,2)
		%alignedWav=shiftVals(goodWaveforms(:,i),goodSpikeShifts(i));
                %wavInterped = spline(1:numBins, alignedWav, 1:splineDt:numBins);
                wavInterped = spline(1:numBins, goodWaveforms(:,i), 1:splineDt:numBins);
                wavInterpedSmooth = fastsmooth(wavInterped, 20, 3, 1);
                %wavInterpedSmooth = wavInterped;
              	%splinedShift=originalBinToSplineBin(goodSpikeShifts(i));
		%wavInterpedSmooth=shiftVals(wavInterpedSmooth,splinedShift);
		goodWaveformsInterp(:,i)=wavInterpedSmooth;
        end

	[dummy, minBins]=min(goodWaveformsInterp,[],1);
	modeMinBin=mode(minBins);
	goodWaveformsInterpTalign=goodWaveformsInterp;
	meanWav=mean(goodWaveformsInterp,2);
	
	for i=1:size(goodWaveformsInterp,2)
		[dummy,minBin]=min(goodWaveformsInterp(:,i));
		alignedWav=shiftVals(goodWaveformsInterp(:,i),modeMinBin-minBin);

		%alignedWav=xcAlignAllPts(goodWaveformsInterp(:,i),meanWav);		

		if(modeMinBin-minBin>50)
			alignedWav=NaN(size(alignedWav));
		end
		%alignedNormedWav=prcNormalize(alignedWav,0,100);
		%alignedNormedWav=minNormalize(alignedWav);
		%goodWaveformsInterpTalign(:,i)=alignedNormedWav;
		goodWaveformsInterpTalign(:,i)=alignedWav;
	end

	plotCheck=0;
	if(plotCheck)
		smallISIs=[diff(cellProp.spikeTimes) 0];
		maxISI=0.05;
		%close all; histogram(smallISIs(smallISIs<maxISI))
		close all; plot(zeros(size(smallISIs(goodIsolatedSpikes))),smallISIs(goodIsolatedSpikes),'ro','Markersize',1)
		title(sprintf('Spikes below %d ms ISI: %d, Total num spikes: %d',maxISI*1000,length(smallISIs(goodIsolatedSpikes)), length(cellProp.spikeTimes)))
		ylim([0 0.1])
		saveas(gcf,'testISIhist.tif')
		[dummy, newMinBins]=min(goodWaveformsInterpTalign,[],1);
		 %close all; histogram(cellProp.spikeMinIndex)
		 close all; histogram(newMinBins)
		 %close all; histogram(goodSpikeShifts)
		saveas(gcf,'testSpikeMinIdxhist.tif')

		close all; plot(goodWaveformsInterpTalign(:,1:20))
		%close all; plot(goodWaveformsInterp)
		meanWav=mean(goodWaveformsInterp,2);
		%meanWav=minNormalize(mean(goodWaveformsInterp,2));
		meanWavAlign=mean(goodWaveformsInterpTalign,2);
		meanWav=troughAlignTo2nd(meanWav,meanWavAlign);

		hold on
		plot(meanWav,'k','LineWidth',7)
		plot(meanWavAlign,'g','LineWidth',7)
		plot((meanWav-meanWavAlign),'r','LineWidth',7)
		
		saveas(gcf,'testAlignedIndividualSpikes.tif')	
		fds
	end
end

function [rsCount,fsCount] = plotSpikeVsNonSpikeLFPDist(cellPropPath,rsCount,fsCount,rsHa,fsHa)

	%load decimated lfp corresponding to this cell's channel
	[decLFP,decFs]=getDecLFPforCell_flux(cellPropPath);

	%zscored and low pass
	%highCut=100;
	%highCut=100;
	%lowpassLFP=filterLFP(decLFP,decFs,0,highCut,2,1);

	%just z-score, no filtering
	lowpassLFP=zscoreLFP(decLFP);

	cellProp=load(cellPropPath);
	spikeTimes=cellProp.spikeTimes(cellProp.goodSpikes);

	if(length(spikeTimes)<25)
		return
	end
	
	%remove spikes that occur within 10 ms of each other
	spikeTimes=removeCloseSpikes(spikeTimes,0.01);

	
	isRS=0;
	isFS=0;
	figNum=0;

	if(cellProp.isInterneuronCell==0)
		rsCount=rsCount+1;
		isRS=1
		figNum=1;
		cellStr='RS';
		if(exist('rsHa'))
			tightHa=rsHa;
		end
	end
	if(cellProp.isInterneuronCell==1)
		fsCount=fsCount+1;
		isFS=1;
		cellStr='FS';
		figNum=2;
		if(exist('fsHa'))
			tightHa=fsHa;
		end
	end

	

	%lookAtIndividualSpikes=1;
	lookAtIndividualSpikes=0;
	if(lookAtIndividualSpikes)
		%dispLFP=decLFP;
		dispLFP=lowpassLFP;
		for j=1:15
			localTime=(spikeTimes(j)-0.005):1/decFs:(spikeTimes(j)+0.01);
			figure(4)
			plot(localTime,dispLFP(round(localTime*decFs)))
			hold on
			plot(spikeTimes(j),dispLFP(round(spikeTimes(j)*decFs)),'r*')
			plot(spikeTimes(j)-1/decFs,dispLFP(round(((spikeTimes(j)-1/decFs)*decFs))),'b*')
		end

		%localTime=(spikeTimes(1)-0.005):1/decFs:(spikeTimes(15)+0.01);
		%plot(localTime,decLFP(round(localTime*decFs)),'k')
		fds

	end

	%minus 1 step (half ms) to avoid rising phase of spike if not smooth!
	%lfpSpikeIdxes=round(spikeTimes*decFs)-1;
	lfpSpikeIdxes=round(spikeTimes*decFs);
	originalRandIdxes=randi(length(lowpassLFP),round(length(spikeTimes)),1);
	%originalRandIdxes=randi(length(lowpassLFP),180000,1);

	duringSpikeLFPvalues=lowpassLFP(lfpSpikeIdxes);

	avoidSpikeTimes=0;
	if(avoidSpikeTimes)
		numRand=0;
		initialSizeFact=2;
		
		while(numRand<length(spikeTimes))

			randIdxes=randi(length(lowpassLFP),round(length(spikeTimes)*initialSizeFact),1);

			minTimeFromSpike=0.1;
			%minTimeFromSpike=0;
			minDistFromNearestSpike=minTimeFromSpike*decFs;
			badIdxes=[];
			tic
			disp('finding indices near spikes')
			for i =1:length(randIdxes)

				if(min(abs(randIdxes(i)-lfpSpikeIdxes))<=minDistFromNearestSpike)	
					badIdxes=[badIdxes i];
				end
			end
			toc

			%size(badIdxes)

			%size(randIdxes)

			randIdxes(badIdxes)=[];

			initialSizeFact=initialSizeFact+1;

			numRand=length(randIdxes);

		end


		randIdxes=randIdxes(1:length(spikeTimes));

	%fds
	%randomLFPvalues=lowpassLFP(randIdxes);
	end

	%0.5ms before spike time to avoid rising phase
	randomLFPvalues=lowpassLFP(originalRandIdxes-1);

	binWidth=(prctile(lowpassLFP,99)-prctile(lowpassLFP,1))/100;
	edges=prctile(lowpassLFP,1):binWidth:prctile(lowpassLFP,99);

	[spikeLFPcounts,edges]=histcounts(duringSpikeLFPvalues,edges);
	[randomLFPcounts,edges]=histcounts(randomLFPvalues,edges);

	%fds

	if(figNum>0)

	figure(figNum)
		if(isRS)
			if(exist('tightHa'))
				axes(tightHa(rsCount))
			end
			%subplot(10,10,rsCount)
		elseif(isFS)
			if(exist('tightHa'))
				axes(tightHa(fsCount))
			end
			%subplot(10,10,fsCount)
		end

		binCenters=edges(1:(end-1))+binWidth/2;
		plot(binCenters,spikeLFPcounts,'r-')
		hold on
		plot(binCenters,randomLFPcounts,'b-')
		%plot(binCenters,randomLFPcounts/max((randomLFPcounts))*max(spikeLFPcounts),'b-')

		%legend('0.5 ms before spike',sprintf('At least %0.1f s from spikes',minTimeFromSpike),'Location','Best')
		%xlabel(sprintf('Lowpass (0-%dHz) LFP z-score',highCut))
		xlabel(sprintf('Raw LFP z-score'))
		ylabel('Count')

	
		set(gca,'fontsize',5)	
		
		numSubs=80;
		axes(tightHa(numSubs))
		plot([0 0],[0 0])
		hold on
		plot([0 0], [0 0])
		set(tightHa(numSubs),'XTickLabel',''); set(tightHa(numSubs),'YTickLabel','')
		if(avoidSpikeTimes)
			%legend('During spike',sprintf('At least %0.1f s from spikes',minTimeFromSpike),'Location','Best')
			lgd=legend(sprintf('0.5 ms before %s spike times',cellStr),sprintf('At least %0.1f s from spikes',minTimeFromSpike),'Location','Best')
		else
			%lgd=legend(sprintf('0.5 ms before %s spike times',cellStr),sprintf('Random timepoints'),'Location','Best')
			lgd=legend(sprintf('Random timepoints'),sprintf('0.5 ms before %s spike times',cellStr),'Location','Best');
		end
		set(gcf,'units','normalized','outerposition',[0 0 1 1])
		drawnow
	end

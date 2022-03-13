function [ax1,ax2,tSgram,fSgram,Sgramlog]=plotRawLFPwithPatternIntensities(ch,singleCycleType,cellLetter)
        disp('Calculating spectrogram............')
	%cellPropFile=load('/nfs/turbo/lsa-ojahmed/processedHumanData/MG49/sessionID-3/cellProperties-MatFiles/06a_cell_properties_MG49.mat');


	chStr=num2str(ch);
	if(ch<10)
		chStr=['0' chStr];
	end

	cellPropFile=sprintf('/nfs/turbo/lsa-ojahmed/processedHumanData/MG49/sessionID-3/cellProperties-MatFiles/%s%s_cell_properties_MG49.mat',chStr,cellLetter);

	lfpFile=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/MG49/sessionID-3/decimatedLFP-MatFiles/Fs2000Hz/concatChan%d/lfpDecConcat.mat',ch);
	
	singleCycleFile=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/MG49/sessionID-3/singleCycleProperties-MatFiles/%s/AlphaSingleCycleBroadbandStats-Ch%s.mat',singleCycleType,chStr);

	lfpData=load(lfpFile);
	lfp=(lfpData.concatenatedLFP);
	Fs=2000;
	
	%maxTime=10000;
	maxTime=length(lfp)/Fs;;
	maxTimeIndex=maxTime*Fs;


	lfp=lfp(1:maxTimeIndex);

	if(strcmp(singleCycleType,'alpha8_2-13Hz'))
		highPass=8.2;
		lowPass=13;
	end

	filteredLFP=filterLFP(lfp,Fs,highPass,lowPass,2,1);

	cellProps=load(cellPropFile);

	spikeTimes=cellProps.spikeTimes;        
	spikeAmps=cellProps.spikeOverallMaxTPAmp;

	

	singleCycleInfo=load(singleCycleFile);

	singleCycleTroughTimes=singleCycleInfo.asymData.fastMidTimes;

	singleCyclePeakToTroughAmps=singleCycleInfo.asymData.fastStartAmp-singleCycleInfo.asymData.fastMidAmp;
	%fds

		
        %% CALCULATE THE SPECTROGRAM
        %spectrumWindow = 30;
        %spectrumWinstep = 30;
         %spectrumWindow = 0.5;
        %spectrumWinstep = 0.1;
         spectrumWindow = 1;
        spectrumWinstep = 1;
         %spectrumWindow = 0.5;
        %spectrumWinstep = 0.5;
        %spectrumWindow = 1;
        %spectrumWinstep = 0.5;
        tapers = [1 1];
        %tapers = [2 3];
        %fpass = [0.5 120];
        fpass = [0.5 20];
        [Sgram fSgram tSgram] = omarSpectrogram(lfp, Fs, spectrumWindow, spectrumWinstep, tapers, fpass, 0);
        Sgramlog = log(Sgram);

    

%         tic
%         window=1;
%         winstep=0.01;
%         lfp=double(lfp);
%         [spectrogram f t]=omarSpectrogram(lfp,Fs,window,winstep);
% 
%         spectrogram=spectrogram';

        figH=figure('units','normalized','outerposition',[0 0 1 1])
        ax1=subplot(4,1,1)
        %plot((1:length(lfp))/Fs,zscoreLFP(lfp))
        reduce_plot((1:length(lfp))/Fs,zscoreLFP(lfp),'b-', (1:length(lfp))/Fs,filteredLFP,'r-')
	ylim([-10 10])
	%hold on
	%plot((1:length(lfp))/Fs,filteredLFP)
	%reduce_plot((1:length(lfp))/Fs,filteredLFP)
	alpha 0.1
	title('time domain intensity')
    xlabel('Time (seconds)');
         legend('Raw LFP', singleCycleType, 'Location','Best')

        ax2=subplot(4,1,3)
            %% PLOT the spectrogram
		omarPcolor(tSgram,fSgram,Sgramlog',figH);
		shading flat
		caxis([-4 15])
		colormap jet;
        %colorbar
        %ylim([0.52 60])
        xlabel('Time (seconds)');
        ylabel('LFP Freqeuncy (Hz)');
        title('Spectrogram')


	ax3=subplot(4,1,4)

	
	isi=[diff(spikeTimes) 0];
	reduce_plot(spikeTimes,1./isi,'ro','MarkerSize',2)
	hold on
	reduce_plot(spikeTimes,1./isi,'r')

	xlabel('Time (s)')
	title('inst. spike rate')

	ax4=subplot(4,1,2)

	reduce_plot(singleCycleTroughTimes,singleCyclePeakToTroughAmps,'k')
	hold on 
	reduce_plot(singleCycleTroughTimes,singleCyclePeakToTroughAmps,'ko','MarkerSize',1)
	xlabel('Time (s)')
	title(sprintf('%s cycle amps',singleCycleType))
	ax1.XAxis.Exponent=0;
	ax2.XAxis.Exponent=0;
	ax3.XAxis.Exponent=0;
	ax4.XAxis.Exponent=0;

	linkaxes([ax1,ax2,ax3,ax4],'x')



	disp('finding cycles with spikes.....')
	singleCycleIndicesWithSpikes=zeros(length(spikeTimes),1);
	for i=1:length(spikeTimes)
		i/length(spikeTimes)
		timeDiffs=abs(spikeTimes(i)-singleCycleTroughTimes);
		[minDummy,minIdx]=min(timeDiffs);
		singleCycleIndicesWithSpikes(i)=minIdx;
	end

	figure
	subplot(1,2,1)
	plot(1./isi,singleCyclePeakToTroughAmps(singleCycleIndicesWithSpikes),'bo','MarkerSize',1)
	subplot(1,2,2)
	randomCycleIndices=randi(length(singleCyclePeakToTroughAmps),size(singleCycleIndicesWithSpikes));
	plot(1./isi,singleCyclePeakToTroughAmps(randomCycleIndices),'bo','MarkerSize',1)



	duringSpikeAlphaAmps=singleCyclePeakToTroughAmps(singleCycleIndicesWithSpikes);
	randomAlphaAmps=singleCyclePeakToTroughAmps(randomCycleIndices);

	
	binWidth=(prctile(singleCyclePeakToTroughAmps,99)-prctile(singleCyclePeakToTroughAmps,1))/100;
        edges=prctile(singleCyclePeakToTroughAmps,1):binWidth:prctile(singleCyclePeakToTroughAmps,99);

	[spikeLFPcounts,edges]=histcounts(duringSpikeAlphaAmps,edges);
        [randomLFPcounts,edges]=histcounts(randomAlphaAmps,edges);	

	figure
	 binCenters=edges(1:(end-1))+binWidth/2;
                plot(binCenters,spikeLFPcounts,'r-')
                hold on
                plot(binCenters,randomLFPcounts,'b-')
	xlabel('Single alpha cycle amplitude')
	legend('Cycles with spikes','Random Cycles')
        %colormap(ax2,fire)
        %startFreq=1;
        %endFreq=100;
        %imagesc(spectrogram,[min(spectrogram(endFreq,:)),max(spectrogram(startFreq,:))])
        %imagesc(spectrogram,[0,1])
        %imagesc(spectrogram)

        
        
        
        
        

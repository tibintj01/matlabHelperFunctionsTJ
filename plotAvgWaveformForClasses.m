function [] = plotAverageWaveformForClass(isInterneuron)

	nexFilePaths=getFilePathsRegex(pwd,'*nex*.mat');
	butterFilePaths=getFilePathsRegex(pwd,'*butter*.mat');
	figure
	timeAxis=(1:48)/30;
	rsWaveform=zeros(48,1);
	fsWaveform=zeros(48,1);
	rsCount=0;
	fsCount=0;
	for cellSzIdx=1:length(isInterneuron)
		cellSzProps=load(nexFilePaths{cellSzIdx});
		goodSpikes=cellSzProps.goodSpikes;

		cellSzProps=load(butterFilePaths{cellSzIdx});

		if(isnan(isInterneuron(cellSzIdx)))
			continue
		end
		if(~isInterneuron(cellSzIdx))
			subplot(2,1,1)	
			%plot(timeAxis,scaledata(cellSzProps.avgGoodWave,-1,1),'r')
			plot(timeAxis,scaledata(nanmean(cellSzProps.waveforms(:,goodSpikes),2),-1,max(nanmean(cellSzProps.waveforms(:,goodSpikes)))),'r')
			hold on
			ylim([-1.2 1.2])
			%rsWaveform=rsWaveform+scaledata(cellSzProps.avgGoodWave,-1,1);
			rsWaveform=rsWaveform+scaledata(nanmean(cellSzProps.waveforms(:,goodSpikes),2),-1,max(nanmean(cellSzProps.waveforms(:,goodSpikes))));
			
			rsCount=rsCount+1;
			title('RS-classified waveforms')
		else
			subplot(2,1,2)	
			%plot(timeAxis,scaledata(cellSzProps.avgGoodWave,-1,1),'b')
			plot(timeAxis,scaledata(nanmean(cellSzProps.waveforms(:,goodSpikes),2),-1,max(nanmean(cellSzProps.waveforms(:,goodSpikes)))),'b')
			hold on
			ylim([-1.2 1.2])
			%fsWaveform=fsWaveform+scaledata(cellSzProps.avgGoodWave,-1,1);
			fsWaveform=fsWaveform+scaledata(nanmean(cellSzProps.waveforms(:,goodSpikes),2),-1,max(nanmean(cellSzProps.waveforms(:,goodSpikes))));
			title('FS-classified waveforms')
			fsCount=fsCount+1;
		end
	end

	
	
	subplot(2,1,1)
	plot(timeAxis,rsWaveform/rsCount,'k','LineWidth',5)

	subplot(2,1,2)
	plot(timeAxis,fsWaveform/fsCount,'k','LineWidth',5)
	xlabel('Time (ms)')
		

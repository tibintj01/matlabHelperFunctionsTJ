function [] = plotAverageWaveformForClass(isInterneuron)

	butterFilePaths=getFilePathsRegex(pwd,'*butter300-Order2.mat');
	figure
	timeAxis=(1:48)/30;
	rsWaveform=zeros(48,1);
	fsWaveform=zeros(48,1);
	rsCount=0;
	fsCount=0;
	for cellSzIdx=1:length(isInterneuron)

		cellSzProps=load(butterFilePaths{cellSzIdx});
		goodSpikes=cellSzProps.goodSpikes;

		if(isnan(isInterneuron(cellSzIdx)))
			continue
		end
		if(isInterneuron(cellSzIdx)==0)
			subplot(2,1,1)	
			%plot(timeAxis,scaledata(cellSzProps.avgGoodWave,-1,1),'r')
			%plot(timeAxis,scaledata(nanmean(cellSzProps.waveforms(:,goodSpikes),2),-1,max(nanmean(cellSzProps.waveforms(:,goodSpikes)))),'r')
			plot(timeAxis,nanmean(cellSzProps.waveforms(:,goodSpikes),2),'r')
			hold on
			%ylim([-1.2 1.2])
			%rsWaveform=rsWaveform+scaledata(cellSzProps.avgGoodWave,-1,1);
			%rsWaveform=rsWaveform+scaledata(nanmean(cellSzProps.waveforms(:,goodSpikes),2),-1,max(nanmean(cellSzProps.waveforms(:,goodSpikes))));
			if(isempty(find(isnan(nanmean(cellSzProps.waveforms(:,goodSpikes),2)))))
				rsWaveform=rsWaveform+nanmean(cellSzProps.waveforms(:,goodSpikes),2);
					
				rsCount=rsCount+1;
			end
			title('RS-classified waveforms')
		elseif(isInterneuron(cellSzIdx)==1)
			subplot(2,1,2)	
			%plot(timeAxis,scaledata(cellSzProps.avgGoodWave,-1,1),'b')
			%plot(timeAxis,scaledata(nanmean(cellSzProps.waveforms(:,goodSpikes),2),-1,max(nanmean(cellSzProps.waveforms(:,goodSpikes)))),'b')
			plot(timeAxis,nanmean(cellSzProps.waveforms(:,goodSpikes),2),'b')
			hold on
			%ylim([-1.2 1.2])
			%fsWaveform=fsWaveform+scaledata(cellSzProps.avgGoodWave,-1,1);
			%fsWaveform=fsWaveform+scaledata(nanmean(cellSzProps.waveforms(:,goodSpikes),2),-1,max(nanmean(cellSzProps.waveforms(:,goodSpikes))));
			fsWaveform=fsWaveform+nanmean(cellSzProps.waveforms(:,goodSpikes),2);
			title('FS-classified waveforms')
			fsCount=fsCount+1;
		end
	end

	
	
	subplot(2,1,1)
	plot(timeAxis,rsWaveform/rsCount,'k','LineWidth',5)

	subplot(2,1,2)
	plot(timeAxis,fsWaveform/fsCount,'k','LineWidth',5)
	xlabel('Time (ms)')
	fds	

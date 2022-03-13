function [] = plotWaveformsByCell(cells,useRaw)

	session_list_human_lfp_spike_relationship_for_students
	sessID=3;
	
	figure
	
	qualThresh=2;
	numCells=length(find(dataInfo(sessID).cellQuality>qualThresh));

	numRows=floor(sqrt(numCells));
	numCols=ceil(sqrt(numCells));
	maxNumWaveformsDisp=500;
	

	plotCount=1;

	marginWidth=0.05;
	%[ha,pos]=tight_subplot(numRows,numCols,0.005,0.005,0.005);
	[ha,pos]=tight_subplot(numRows,numCols,marginWidth,marginWidth,marginWidth);
	for i=1:length(cells)
		cell=cells{i};
		if(cell.qualityRating>qualThresh && plotCount<=numRows*numCols)
			%subplot(numRows,numCols,plotCount)
			axes(ha(plotCount))
			
			%default is NEX waveForms
			waveforms=cell.waveForms;
			numWavs=min(maxNumWaveformsDisp,size(waveforms,1));

			if(useRaw==2)
				waveformsRaw=cell.rawSpikeWaveforms(1:numWavs,:);
                                startIdx=cell.spikeIdxInRaw-9;
                                finishIdx=cell.spikeIdxInRaw+38;

                                %waveforms=waveformsRaw(:,startIdx:finishIdx);

				filtWaveforms=getFilteredWaveforms(waveformsRaw,cell.Fs);
				waveforms=filtWaveforms(:,startIdx:finishIdx);
			elseif(useRaw==1)	
				waveformsRaw=cell.rawSpikeWaveforms(1:numWavs,:);
				startIdx=cell.spikeIdxInRaw-9;
				finishIdx=cell.spikeIdxInRaw+38;
				
				waveforms=waveformsRaw(:,startIdx:finishIdx);
				medians=median(waveforms,2);
				waveforms=waveforms-medians;

			end

			wavesDisp=waveforms(1:numWavs,:);

			
			plotCount=plotCount+1;
			plot(wavesDisp')
			title(i)
			%title(dataInfo(sessID).cellIsInterneuron(i))
			%ylim([-0.25, 0.25])
			xlim([1 48])
			set(gca,'FontSize',3)
			drawnow
		end
	end
	set(gcf,'units','normalized','outerposition',[0 0 1 1])

	addStr='';
	if(useRaw==1)
		addStr='raw';
	end

	if(useRaw==2)
		addStr='butterworth600To3000Hz';
	end
	
	saveas(gcf,sprintf('cellWaveformsPlot_%s.tif',addStr))


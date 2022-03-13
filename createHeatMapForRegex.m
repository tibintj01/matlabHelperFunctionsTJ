function [] = createHeatMapForRegex(dirName,regexFileName)
	filePaths=getRegexFilePaths(dirName,regexFileName);

	[lowFreq,highFreq]=getFreqBand(regexFileName);

	testData=load(filePaths{1});

	periSpikeHeatMap=zeros(length(filePaths),length(testData.periTime));

	for i=1:length(filePaths)
		data=load(filePaths{i});
		periSpikeHeatMap(i,:)=data.periSpikeLFP;
	end
	[minHeat,minHeatIdx]=min(periSpikeHeatMap,[],2);
	[meanSorted,meanOrder]=sort(minHeatIdx);
	periSpikeHeatMap=periSpikeHeatMap(meanOrder,:);

	disp('plotting heat map....')
	figure
	colormap(jet)
	imagesc(periSpikeHeatMap)
	colorbar

	disp('setting x labels.....')
	setXLabelsImage(periSpikeHeatMap,[min(testData.periTime),0,max(testData.periTime)]);

	ylabel('Cell No.')
	xlabel('Time after spike (s)')

	ylabel(colorbar,'LFP')	
	if(contains(regexFileName,'RS'))
		title(sprintf('RS cells; avg. peri-spike LFPs, filtered %.2f to %2.f Hz',lowFreq,highFreq))
	elseif(contains(regexFileName,'FS'))
		title(sprintf('FS cells; avg. peri-spike LFPs, filtered %.2f to %2.f Hz',lowFreq,highFreq))
		%title('FS cells')

	elseif(contains(regexFileName,'GS'))
		title(sprintf('Gray cells; peri-spike LFPs, filtered %2.f to %2.f Hz',lowFreq,highFreq))
		%title('Gray cells')

	end
	saveName=strrep(regexFileName(1:(end-4)),'*','-');	
	saveas(gcf,fullfile(dirName,[saveName '.tif']))

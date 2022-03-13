
		saveDir='/nfs/turbo/lsa-ojahmed/periSpikeFilteredLFPs_MG49';
		if(~isdir(saveDir))
			mkdir(saveDir)
		end

       parpool(7) 
        %for i=1:length(freqBands)
        parfor i=1:7
		cellPropDir='/nfs/turbo/lsa-ojahmed/classifiedMG49-sleepWake'
		saveDir='/nfs/turbo/lsa-ojahmed/periSpikeFilteredLFPs_MG49';
		freqBands=freqBandNames();
		cellPropFileNames=getRegexFilePaths(cellPropDir,'*cell_prop*mat');
		maxNumSpikes=1000000;
		count=0;
                for j=1:length(cellPropFileNames)
                        %saveAllPeriSpikeLFPforCell(cellPropFileNames{j},freqBands{i},maxNumSpikes)
                        saveAllPeriSpikeLFPforCell(cellPropFileNames{j},freqBands{i},maxNumSpikes,saveDir)
                end
        end

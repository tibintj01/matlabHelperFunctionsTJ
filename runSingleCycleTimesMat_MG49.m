
		%saveDir='/nfs/turbo/lsa-ojahmed/mg49_singleCycleBroadBandStats_Par';
		saveDir='/nfs/turbo/lsa-ojahmed/processedHumanData/MG49/sessionID-3/singleCycleProperties-MatFiles/alpha6-40Hz';
		saveDir='/nfs/turbo/lsa-ojahmed/processedHumanData/MG49/sessionID-3/singleCycleProperties-MatFiles/alpha8_2-13Hz';
		if(~isdir(saveDir))
			mkdir(saveDir)
		end

       parpool(19) 
        %for i=1:length(freqBands)
        %for i=1:96
        parfor i=1:96
		%saveDir='/nfs/turbo/lsa-ojahmed/mg49_singleCycleBroadBandStats_Par';
		saveDir='/nfs/turbo/lsa-ojahmed/processedHumanData/MG49/sessionID-3/singleCycleProperties-MatFiles/alpha8_2-13Hz';
		%freqBands=freqBandNames();
		freqBands={'Alpha'};
                for j=1:length(freqBands)
                        %saveAllPeriSpikeLFPforCell(cellPropFileNames{j},freqBands{i},maxNumSpikes)
                      	try
			  saveAllLFPsingleCyclesForFreqBand(i,freqBands{j},saveDir)
			catch
			   disp(sprintf('error in getting single cycles, skipping ch %d for freqBand %s',i,freqBands{j}))
			end
                end
        end

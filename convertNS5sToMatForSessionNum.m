function [] = convertNS5toMatForSessionNum(currSes)

	 %load session meta data
	 session_list_human_lfp_spike_relationship_for_students
	
	 %set to 1 -> will not overwrite existing mat files
	 checkBeforeConverting=1;

	 %inputDataDir='C:/Users/tibin/Desktop/FOR_TIBIN';
	 %outputDir=sprintf('session%dMat',currSes)
	
	 %inputDataDir='/nfs/turbo/lsa-ojahmed/tibin/FOR_TIBIN_otherPts/otherPtData_12_7_17';
	 inputDataDir='/nfs/turbo/lsa-ojahmed/tibin/FOR_TIBIN_otherPts';
	 outputDir='/nfs/turbo/lsa-ojahmed/tibin/processedHumanData';


	%% NS5 PROCESSING: Make sure that the relevant NS5 files have been converted to, for example, lfp6_original.mat and lfp6_dec.mat
	for currNS5Ind = 1:dataInfo(currSes).numNS5Files % LOOP AROUND ALL NS5 FILES
	    currNS5Filename = dataInfo(currSes).ns5Filenames{currNS5Ind};
	    humanConvertNS5toMat(checkBeforeConverting, dataInfo(currSes).subject, currNS5Filename, ...
		dataInfo(currSes).ns5VisChannels, dataInfo(currSes).ns5SeqChannels,inputDataDir,outputDir);
	    if ~isempty(dataInfo(currSes).ns5TrigVisChannel)
		humanConvertNS5toMat(checkBeforeConverting, dataInfo(currSes).subject, currNS5Filename, ...
		    dataInfo(currSes).ns5TrigVisChannel, dataInfo(currSes).ns5TrigSeqChannel,inputDataDir,outputDir);
	    end
	end


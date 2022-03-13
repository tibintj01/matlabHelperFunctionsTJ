function [] = createConcatenatedLFPforSessionNum(sessionIdx,lfpMatDir)
	session_list_human_lfp_spike_relationship_for_students

	currSes=sessionIdx;

	%adapted from Omar's human_GRAND_preprocessing_script.m  
        for i=dataInfo(currSes).ns5VisChannels
		currOutputFilename = sprintf('%s_ch%d', dataInfo(currSes).descriptiveFilename, i);
		omarCombineLFPChanIntoNS5_Tibin(dataInfo(currSes).subject, i, dataInfo(currSes).numNS5Files, dataInfo(currSes).ns5Filenames, currOutputFilename, dataInfo(currSes).descriptiveFilename, 0,lfpMatDir, [lfpMatDir '-concatenated']);
	end


function [concatenatedLFP] = getConcatenateLFPForChan(sessNum,dataDir,chan)

	%dataDir should be directly above recording directories
	 session_list_human_lfp_spike_relationship_for_students

	subject=dataInfo(sessNum).subject;
	 recordingDirs=dataInfo(sessNum).ns5Filenames;

		 concatenatedLFP=[];
		 for fileNum=1:length(recordingDirs)
			recordingDir=recordingDirs{fileNum};
			 filePath=fullfile(dataDir,subject,recordingDir,sprintf('lfp%d_original.mat',chan));
		 	if(exist(filePath))
		 	 	lfpData=matfile(filePath);
			 else
			 	error('channel not found....')
			 end
			lfp=round(double(lfpData.lfp)*0.249*10000)/10000;%convert digital to uV  - see humanProcesLFP... script
			concatenatedLFP=[concatenatedLFP lfp(:)'];
		 end
		


	

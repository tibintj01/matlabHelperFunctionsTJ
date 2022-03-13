parpool(20)
parfor i=1:96
	%sessionID=6
	%sessionID=9 %MG67 includes tasks...... sorted cell list not updated?

	%for sessionID=2:2
	%for sessionID=13:13
	for sessionID=3:3
		try
			if(sessionID==-1)
				convertNS5sToMatForSessionNum(sessionID)
				createConcatenatedLFPforSessionNum(sessionID,lfpMatDir)
			end
			saveDecLFPAndCellPropMatFromConcatNS5Ch(sessionID,i)
		catch ME
			disp(sprintf('error for ch %d: ',i))
			disp(ME.message)
		end
	end
end

%{
session_list_human_lfp_spike_relationship_for_students
ptDir=dataInfo(sessionID).subject;
cellPropSaveDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/%s/sessionID-%d/cellProperties-MatFiles',ptDir,sessionID);
classifyFromCellProp_TibinPtSleepWake(cellPropSaveDir)
%}

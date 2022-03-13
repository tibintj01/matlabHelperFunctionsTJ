function [] = runCellClassificationForSessionID(sessionID)
	%sessionID=6
	session_list_human_lfp_spike_relationship_for_students
	ptDir=dataInfo(sessionID).subject;
	cellPropSaveDir=sprintf('/nfs/turbo/lsa-ojahmed/tibin/processedHumanData/%s/sessionID-%d/cellProperties-MatFiles',ptDir,sessionID);
	%classifyFromCellProp_TibinPtSleepWake(cellPropSaveDir,sessionID)
	classifyFromCellProp_Tibin(cellPropSaveDir,sessionID)

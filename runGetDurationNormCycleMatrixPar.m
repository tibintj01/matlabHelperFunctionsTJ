
parpool(19)
parfor ch=1:96 
	try
		saveDir='/nfs/turbo/lsa-ojahmed/processedHumanData/MG49/sessionID-3/durationNormalizedSingleCycles-MatFiles/broadbandAlpha';
		ptDir='MG49';
		sessionNum=3;
		singleCycleType='broadbandAlpha';
		singleStatsDirToDurationNormedCycleMatrixForCh(ptDir,sessionNum,ch,singleCycleType,saveDir)
	catch
		disp('error...skipping ch')
	end

end

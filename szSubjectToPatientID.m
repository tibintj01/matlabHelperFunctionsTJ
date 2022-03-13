function [ptID]=szSubjectToPatientID(szSubjectName)
	
	ptIDs={'MG49','MG63','BW9','RIHE1'};

	for i=1:length(ptIDs)
		if(contains(szSubjectName,ptIDs{i}))
			ptID=ptIDs{i};
			return
		end
	end

	ptID=NaN;

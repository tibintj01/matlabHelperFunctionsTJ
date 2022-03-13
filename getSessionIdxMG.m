function [sessIdx]=getSeizureStartTimes(cellPropFileName)


	sessIdx='error';

	if(length(findstr(cellPropFileName,'MG29'))>0)
		sessIdx=1;
	end	
	if(length(findstr(cellPropFileName,'MG49'))>0)
		sessIdx=3;
	end	




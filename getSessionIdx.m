function [sessIdx]=getSeizureStartTimes(cellPropFileName)


	sessIdx='error';

	if(length(findstr(cellPropFileName,'MG49_seizure36'))>0)
		sessIdx=1;
	end	
	if(length(findstr(cellPropFileName,'MG49-seizure43'))>0)
		sessIdx=2;
	end	
	if(length(findstr(cellPropFileName,'MG49-seizure45'))>0)
		sessIdx=3;
	end	
	if(length(findstr(cellPropFileName,'MG63_seizure1-4'))>0)
		sessIdx=4;
	end	
	if(length(findstr(cellPropFileName,'BW9_seizure1-3'))>0)
		sessIdx=5;
	end	
	if(length(findstr(cellPropFileName,'RIHE1'))>0)
		sessIdx=6;
	end	




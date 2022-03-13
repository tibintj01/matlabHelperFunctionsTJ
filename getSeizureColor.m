function [szColor]=getSeizureStartTimes(cellPropFileName)


	szColor='error';

	if(length(findstr(cellPropFileName,'MG49_seizure36'))>0)
		szColor='r';
	end	
	if(length(findstr(cellPropFileName,'MG49-seizure43'))>0)
		szColor='g';
	end	
	if(length(findstr(cellPropFileName,'MG49-seizure45'))>0)
		szColor='b';
	end	
	if(length(findstr(cellPropFileName,'MG63_seizure1-4'))>0)
		szColor='k';
	end	
	if(length(findstr(cellPropFileName,'BW9_seizure1-3'))>0)
		szColor='y';
	end	
	if(length(findstr(cellPropFileName,'RIHE1'))>0)
		szColor='m';
	end	




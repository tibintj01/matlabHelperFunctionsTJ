function [szStartTimes]=getSeizureStartTimes(cellPropFileName)

	session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ

	szIdx=0;

	if(length(findstr(cellPropFileName,'MG49_seizure36'))>0)
		szIdx=1;
	end	
	if(length(findstr(cellPropFileName,'MG49-seizure43'))>0)
		szIdx=2;
	end	
	if(length(findstr(cellPropFileName,'MG49-seizure45'))>0)
		szIdx=3;
	end	
	if(length(findstr(cellPropFileName,'MG63_seizure1-4'))>0)
		szIdx=4;
	end	
	if(length(findstr(cellPropFileName,'BW9_seizure1-3'))>0)
		szIdx=5;
	end	
	if(length(findstr(cellPropFileName,'RIHE1'))>0)
		szIdx=6;
	end	



	szStartTimes=dataInfo(szIdx).szStartTime;

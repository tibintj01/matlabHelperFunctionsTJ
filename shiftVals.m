function [shiftedVals] = shiftVals(vals,i)
	vals=vals(:)'; %convert to row vector
	circShifted=circshift(vals,[0,i]);

	if(i>0)
		%if shifting right, replace emptys w/ first value
		circShifted(1:i)=vals(1); 
	elseif(i<0)
		%if shifting left, replace emptys w/ last value
		circShifted((end+i+1):end)=vals(end);
	end
	
	shiftedVals=circShifted;

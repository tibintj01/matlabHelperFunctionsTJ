function [interpedVals,newFs]=getInterped(vals,interpFact,Fs)
	%Tibin John, 9/6/17
	%Description: gets interpolated values series
	%Input: values, interpolation factor
	%Output: interpolated vector

	x=1:length(vals);
	newStep=1/interpFact;
	xQueries=1:newStep:length(vals);
	
	newFs=Fs*interpFact;
	
	interpedVals=interp1(x,vals,xQueries,'spline');

	%interpedVals=zeros(size(xQueries));
	%count=1;
	%for xq=xQueries
	%	interpedVals(count)=interp1(x,vals,xq,'spline');
	%	count=count+1;
	%end

	

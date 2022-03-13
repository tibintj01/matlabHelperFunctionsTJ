function [szIdx] = getSzIdx(propPath)
        %assumes file saved as *seizure%d.mat
	if(isempty(findstr(propPath,'butter')))
        	szIdx=str2num(propPath(end-8));
	else
        	%szIdx=str2num(propPath(end-14));
        	szIdx=str2num(propPath(end-21));
	end

function [strIdx] = getIdxInStrCellArray(strArray,desiredStr)
	[isMem,strIdx]=ismember(desiredStr,strArray);
	if(strIdx==0)
		disp('Pt id name not found.')
	end	

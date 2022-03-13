function [valueIdxes]=getIdxOf(desiredValues,valueArray)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%       Function to get idxes of unique values in array (works for floats)
%Result:
%       Returns array of idxes corresponding to value array provided
%Author: 
%       Tibin John, June 2017, tibin.john015@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%corrects if parameters are given in wrong order
if(length(desiredValues)>length(valueArray))
	temp=valueArray(:);
	valueArray=desiredValues(:);
	desiredValues=temp(:);
end

for i=1:length(desiredValues)
        [dummyVals,valueIdxes(i)]=min(abs(valueArray-desiredValues(i)));
end

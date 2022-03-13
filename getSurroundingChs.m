function [neighborChs] = getSurroundingChs(ch,ptDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%
%
%Input:
%
%
%Output:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-06-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[arrayMap electrodeXY electrodeImp] = neuroportArrayData(ptDir);
arrayMap=fliplr(arrayMap)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')
col=electrodeXY(ch,1);
row=(10-electrodeXY(ch,2))+1;
%row=(electrodeXY(ch,2));

neighborChs=NaN(9,1);

count=1;
for horiz=-1:1
   for vert=-1:1
	if(row+vert<=size(arrayMap,1) && col+horiz<=size(arrayMap,2) && row+vert>=1 && col+horiz>=1)
        	neighborChs(count)=arrayMap(row+vert,col+horiz);
		count=count+1;
	end
   end
end
neighborChs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

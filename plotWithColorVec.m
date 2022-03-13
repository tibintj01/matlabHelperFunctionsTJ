function [] = plotWithColorVec(x,y,c,colorMapName,lineWidth,zLevel)
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
%Created on 2018-09-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing plotWithColorVec fxn output.........')

z=zeros(size(x));
if(exist('zLevel'))
	z=z+zLevel;
end

if(~exist('lineWidth'))
	h=surface([x;x],[y;y],[z;z],[c;c],'facecol','no','edgecol','flat','linew',1);
else
	h=surface([x;x],[y;y],[z;z],[c;c],'facecol','no','edgecol','flat','linew',lineWidth);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%colormap jet

colormap(gca,eval(colorMapName))


%if(strcmp(colorMapName,'copper'))
%	currCaxis=caxis;
%	caxis=[currCaxis(1) currCaxis(2)*0.9];
%end

colorbar

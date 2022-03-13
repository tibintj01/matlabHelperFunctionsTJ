function [pointCloud]=getPointCloudFromWaveforms(waveforms,U)
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
%Created on 2018-09-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numFeatures=4;
pointCloud=NaN(size(waveforms,2),numFeatures);
for waveNum=1:size(waveforms,2)
	currWaveform=waveforms(:,waveNum);
	[pc1,pc2]=projectWaveformOntoU(currWaveform,U);
	%[currWaveformTPwidth,currWaveformTPamp]=getWidthAmpOfWaveform(currWaveform);
	[currWaveformTPwidth,currWaveformTPamp]=getTPwidthAmpOfWaveform(currWaveform);

	pointCloud(waveNum,:)=[pc1 pc2 currWaveformTPwidth currWaveformTPamp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing getPointCloudFromWaveforms fxn output.........')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [waveformHistMatrix,tBins,vBins] = getSpikeHistHeapMap(waveformsInColumns,minV,maxV)
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
%Created on 2018-09-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('minV','var'))
	minV=-200;
	maxV=100;
end
numVbins=200;
vEdges=linspace(minV,maxV,numVbins+1);
vBins=edgesToBins(vEdges);

numTbins=size(waveformsInColumns,1)
tBins=1:numTbins;

waveformHistMatrix=zeros(numVbins,numTbins);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing getSpikeHistHeapMap fxn output.........')

for i=1:numTbins
	tBinVhist=histcounts(waveformsInColumns(i,:),vEdges);
	waveformHistMatrix(:,i)=tBinVhist;	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


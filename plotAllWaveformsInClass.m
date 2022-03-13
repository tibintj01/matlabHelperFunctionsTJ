function [outputStruct] = plotAllWaveformsInClass(RS0FS1Int2,cellPropDir,cellPropIDStr)
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
%Created on 2018-08-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs=30000;
%Fs=32000.0
ephysFigureColors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract input variables from struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RS0FS1Int2=inputStruct.RS0FS1Int2;

%try
%	cellPropDir=inputStruct.cellPropDir;
%catch


if(RS0FS1Int2==0)
	cellColor=RSColor;
elseif(RS0FS1Int2==1)
	cellColor=FS1Color;
else
	cellColor=FS1ColorLight;
end

if(~exist('cellPropDir','var'))
	cellPropDir='/nfs/turbo/lsa-ojahmed/tibin/tibinCellSzClassificationFiles';
end

if(~exist('cellPropIDStr','var'))
	allCellPropFiles=getRegexFilePaths(cellPropDir,sprintf('*cell_prop*mat'));
else
	allCellPropFiles=getRegexFilePaths(cellPropDir,sprintf('*cell_prop*%s*mat',cellPropIDStr));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing plotAllWaveformsInClass fxn output.........')

close all
for i=1:length(allCellPropFiles)
	currCellProps=load(allCellPropFiles{i});
	if(currCellProps.isInterneuronCell~=RS0FS1Int2)
		continue
	end

	currWaveform=currCellProps.avgGoodWave;
	timeAxis=(1:length(currWaveform))/Fs*1000;
	currWaveformNorm=currWaveform/(abs(min(currWaveform)));
	%plot(timeAxis,currWaveformNorm,'Color',cellColor)
	if(i==2)
	plot(timeAxis,currWaveform,'Color',cellColor)
	hold on	
	i/length(allCellPropFiles)	

	title(removeUnderscores(allCellPropFiles{i}))	
	allCellPropFiles{i}
	saveas(gcf,'checkAmp.tif')
	fds
	end
end

xlabel('Time (ms)')
ylabel('A.U.')
xlim([-Inf Inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store output variables in struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


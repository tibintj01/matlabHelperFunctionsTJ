function [cellWaveformTiming] = getCellWaveformTimingMatrixForPt(ptID)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Tibin John, 8/13/17
	%Description: Accumulates NeuroExplorer spiking data across channels of electrode array for a given patient
	%Input: patient ID string
	%Output:one cell array of structs for each cell, containing timing, waveforms, and other info
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%load meta data for all available cells 

	session_list_human_lfp_spike_relationship_for_students
	%Loaded variables:
	% dataInfo (1x18 struct)
	%

	%get patient idx from ID string
	%multiple entries per patient (various tasks/concatenations) - correct for this!
	if(~exist('ptID'))
		ptID='MG49'
	end
	ptIDs={dataInfo.subject};
	sessionNum=getIdxInStrCellArray(ptIDs, ptID)


	%get channel corresponding to each cell detected
	cellChannels=dataInfo(sessionNum).cellChannel;
	numCells=length(cellChannels);

	%populate cell array containing spike times, waveforms, and chanNum for ea. cell
	cellWaveformTiming=cell(numCells,1);
	for cellNum=1:length(cellChannels)
		chanNum=cellChannels(cellNum);
		nexFileName=getNexFileName(chanNum);
		cellIdxInNex=getCellIdxWithinChannel(cellNum,cellChannels);
		nexData=readNexFile(nexFileName);
	 	cellWaveFormTiming{cellNum}=struct();

		cellWaveformTiming{cellNum}.spikeTimes=nexData.neurons{cellIdxInNex}.timestamps;
		%row numbers represent spike number, columns represent increasing time
		cellWaveformTiming{cellNum}.waveForms=nexData.waves{cellIdxInNex}.waveforms';

		cellWaveformTiming{cellNum}.chanNum=chanNum;
		cellWaveformTiming{cellNum}.qualityRating=dataInfo(sessionNum).cellQuality(cellNum);
		cellWaveformTiming{cellNum}.cellIdxInChannel=cellIdxInNex;
		cellWaveformTiming{cellNum}.Fs=nexData.freq;
		cellWaveformTiming{cellNum}.tbeg=nexData.tbeg;
		cellWaveformTiming{cellNum}.tend=nexData.tend;
	end





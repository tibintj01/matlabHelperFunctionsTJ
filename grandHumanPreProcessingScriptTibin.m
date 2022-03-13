%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description: grandHumanPreProcessingScriptTibin
%
%
%Preconditions:
%	NS5_NEV folders in input dir are populated with ns5 and nev files according to session list from experiment
%	
%Effects:
%	Populates processedHumanDirectory according to directory structure processedHumanDataDirectoryStructure.tif
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-06-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session_list_human_lfp_spike_relationship_for_students
%sessionIndices=[ 1 6 9]; %MG29 MG63 MG67
%sessionIndices=[ 6 9 1]; %MG29 MG63 MG67
sessionIndices=[ 1 2 ]; %MG29 MG63 MG67

for i=1:length(sessionIndices)

	sessionIdx=sessionIndices(i);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Task #1
	%Description:
	% save ns5's as mats according to session list
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('computing output.........')
	lfpMatDir='/nfs/turbo/lsa-ojahmed/tibin/processedHumanData'

	%best on single 500GB node
	if(sessionIdx==-1)
		convertNS5sToMatForSessionNum(sessionIdx)
		createConcatenatedLFPforSessionNum(sessionIdx,lfpMatDir)
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Task #2
	%Description:
	%    save concatenated lfp files in directory corresponding to directory structure
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('computing output.........')
	

	saveDecLFPAndCellPropMatFromConcatNS5AllCh(sessionIdx)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Task #3
	%Description:
	%    get cell properties and add classification fields to cell prop files corresponding to directory structure
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('computing output.........')
	runCellClassificationForSessionID(sessionIdx)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Task #4
	%Description:
	%    generate and save spectrogram segments (copper color map?) 
	%	 and per-frequency analyses in directory correpsonding to directory structure
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('computing output.........')


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Task #5
	%Description:
	%    generate and save LFP-LFP and LFP-spike kappa/PPC vs space plots    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('computing output.........')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Task #6
	%Description:
	%    show spike cross correlations confirming E, I classification
	%    	generate and save spatial extent of spike train correlations for E,I; sleep and wake states
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('computing output.........')


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Task #7
	%Description:
	%    for time segments of interest, generate locking/timing movies of local lfp and spikes
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('computing output.........')


end

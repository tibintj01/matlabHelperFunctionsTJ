function [nexFileName] =getNexFileName(chanNum)
	%returns file path for channel of current patient MG49 on Tibin's computer
	%%%update to generalize for patient and NeuroExplorer file directory%%
	baseStr='C:/Users/tibin/Desktop/FOR_TIBIN/MG49/sorted/20110615-094530-023_024_025_026_030_031_032_plus_anesthesia/chan';
	nexFileName=sprintf([baseStr '%d/20110615-094530-023_024_025_026_030_031_032_plus_anesthesia_ch%d.nex'],chanNum,chanNum);	

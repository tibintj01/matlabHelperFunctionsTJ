function [cellTypeStr] =getCellTypeStr(sessionIdx,ch,indWithinCh,szNum)
	%Description:
	%Given session index, channel, index within channel, and seizure number (1-3),
	%returns property structure now with fields isInterneuronCell, isRSprobGMMcell,mahalanobisDistancesCell
	%isInterneuronCell contains 0 if gaussian mixture model posterior probability of being in RS cluster is less than 0.05, 1 if greater than 0.95, and 2 if between 0.05 and 0.95
	%isRSprobGMMcell is gaussian mixture model posterior probability of being in RS cluster
	%mahalanobisDistancesCell is a 2-element array; first number is mahalanobis distance from RS cluster, second is mahalanobis distance from FS cluster
	%props structure also contains all waveforms butterworth 300-3000Hz filtered (2nd order) in the waveforms field (index with goodSpikes to get just single, aligned spikes used to calculate feature vectors) to visually confirm classification
	%also returns descriptive file name containing this property struct as optional second output

	%Tibin John, Nov. 2017


	%initialize input directory containing mat files with class information
	%**CHANGE this to fit the directory containing the *butter300-Order2.mat files**
	classificationPropsDir='/nfs/turbo/lsa-ojahmed/tibin/tibinCellSzClassificationFiles';

	session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ

	%find file name corresponding to input while looping through session list	

	if(sessionIdx==2)
                subjectName='MG49-seizure43';
        else
                subjectName=dataInfo(sessionIdx).subject;
        end
        regexStr=sprintf('%d%s_cell_properties_%s*seizure%d*butter300-Order2.mat',ch,num2letter(indWithinCh),subjectName,szNum);

	filePaths=getRegexFilePaths(classificationPropsDir,regexStr);

	%if file does not exist, this cell was not included in analysis (e.g. quality < 3)
	if(isempty(filePaths))
		cellProps=NaN;
		filePath=NaN;
		cellTypeStr='Unclassified';
		return
	end
	filePath=filePaths{1};
	
	cellProps=load(filePath);	

	if(cellProps.isInterneuronCell==0)
		cellTypeStr='RS';
	elseif(cellProps.isInterneuronCell==1)
		cellTypeStr='FS';
	elseif(cellProps.isInterneuronCell==2)
		cellTypeStr='Intermediate';
	end


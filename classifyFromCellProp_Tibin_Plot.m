function [] = classifyFromCellProp_Tibin_Plot(cellPropDir,ax1,ax2) 
%Description: This code plots combinations of features from cell property files %created by getCellWaveformPropsOmar_Tibin.m
%Tibin John

Fs=30000;
reclassify=0;
session_list_human_lfp_spike_relationship_for_students
%ptDir=dataInfo(sessionID).subject;

%586 cells - different segmentations....
%by seizure, by patient, all together
%2D pairs of features for each

%make table of 6 features along with patient ID and seizure ID (8 total features), ea. row is a cell
%splitapply plot and plot3 based on ptID and szID groupings


%just 40 configurations so just read each from file each time
%All cell:
%586 row (cell), 6 column (feature) array  
%By Patient:
%e.g. 96 row (avg'd over seizures?-weighted by # of spikes!), 6 column array
%By Seizure:
%e.g. 96 row, 6 column array


%segment by directory and file names.... bottom up or top down?

%nexFilePaths=getFilePathsRegex(pwd,'*butter300-Order2.mat');
nexFilePaths=getFilePathsRegex(cellPropDir,'*cell_prop*.mat');

%f1=getFeatureVectorFromNexPropFile(nexFilePaths{1});
[f1,varNames]=getFeatureVectorFromButterOrder2(nexFilePaths{1});

%cellFeatureTable=cell2table(f1,'VariableNames',{'PtID','SzNum','nexTPwidth','rawTPwidth','butter300TPwidth','nexP50','nexP25','nexT50','rawT15','butter300T15','nexT15','butter300P25','butter300T50'});
cellFeatureTable=cell2table(f1,'VariableNames',varNames);

disp('extracting features from cells....')

for fileNum=2:length(nexFilePaths)
	%tic
	%f=getFeatureVectorFromNexPropFile(nexFilePaths{fileNum});
	[f,dummy]=getFeatureVectorFromButterOrder2(nexFilePaths{fileNum});
	cellFeatureTable=[cellFeatureTable;f];
	%toc
	%if(fileNum==25)
	%	fds
	%end
end
%cellFeatureTable

%k means based on all cell-seizure pairs in big 3 features
%isInterneuron=kmeans(cellFeatureTable{:,[5,6,10]},2)-1;
%isInterneuron=kmeans(cellFeatureTable{:,[4,5,8]},2)-1;
%isInterneuron=kmeans(cellFeatureTable{:,[12,24,22]},2)-1;
%isInterneuron=kmeans(cellFeatureTable{:,[12,24]},2,'Replicates',1000)-1;
sigLev=0.05;

%175 is noise, not cell - with Omar Nov. 29
%cellFeatureTable(175,3:end)=array2table(NaN(size(cellFeatureTable(175,3:end))));
%file name corresponding to non-cell:36b_cell_properties_RIHE1_seizure2_butter300-Order2.mat
%nexFilePaths{175}
%badCellProp=load(nexFilePaths{175})
%mean spike T15 width = 22.877 bins, 0.76 ms
%fds

%CONFIRM 12 and 24 are right features
%[isInterneuron,isRSprobGMM,mahalanobisDistances]=getGMMclusterMemberships(cellFeatureTable{:,[12,24]},2,sigLev,ax1);
[isInterneuron,isRSprobGMM,mahalanobisDistances]=getGMMclusterMemberships(cellFeatureTable{:,[12,24]},2,sigLev);
%[isInterneuron,isRSprobGMM,mahalanobisDistances]=getGMMclusterMemberships(cellFeatureTable{:,[16,28]},2,sigLev);
%fds
for f=1:length(nexFilePaths)
	isInterneuronCell=isInterneuron(f);
	isRSprobGMMcell=isRSprobGMM(f);
	mahalanobisDistancesCell=mahalanobisDistances(f,:);
	nexFilePaths{f}
	isInterneuronCell
	if(reclassify)
		save(nexFilePaths{f},'isInterneuronCell','-append');
		save(nexFilePaths{f},'isRSprobGMMcell','-append');
		save(nexFilePaths{f},'mahalanobisDistancesCell','-append');
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

replot=1;
plotGroups=0;

if(replot)
	groupName='All_Patients_and_seizures_spike_features';
	%groupName=sprintf('%s-SessionID-%d',ptDir,sessionID);
		%plotFeatures2D(cellFeatureTable,'spikeT15WidthButter300Order2','spikeTPWidthButter300Order2',groupName,isInterneuron)
		%plotFeatures2D(cellFeatureTable,'spikeT15WidthButter300Order2','spikeTPWidthButter300Order2',groupName,isInterneuron,[],Fs,ax2)
		plotFeatures2D(cellFeatureTable,'spikeT15WidthButter300Order2','spikeTPWidthButter300Order2',groupName,isInterneuron,[],Fs,[])
		%plotFeatures3D(cellFeatureTable,'spikeT15WidthButter300Order2','spikeTPWidthButter300Order2','spikeP75WidthButter300Order2',groupName,isInterneuron)
	%plotAllCombinations(cellFeatureTable,[],groupName)
	%plotFeatures3D(cellFeatureTable,'butter300T15','butter300TPwidth','nexP50',groupName,isInterneuron)
	%plotFeatures3D(cellFeatureTable,'butter300T15','nexP25','butter300P25',groupName)
	%plotFeatures3D(cellFeatureTable,'butter300T15','butter300TPwidth','butter300P25',groupName)
	%plotFeatures3D(cellFeatureTable,'butter300T15','nexT50','nexP25',groupName)
	%plotFeatures3D(cellFeatureTable,'butter300T15','butter300TPwidth','butter300TFirstPwidth',groupName,isInterneuron)
	%fds

	%close all

	if(plotGroups)

	[ptGroupIdxes,ptGroupNames]=findgroups(cellFeatureTable.PtID);

	ptGroupIdxes=ptGroupIdxes(~isnan(ptGroupIdxes));

	for ptIdx=1:length(unique(ptGroupIdxes))
		groupName=sprintf('Patient-%s',ptGroupNames{ptIdx})
		%plotFeatures3D(cellFeatureTable(ptGroupIdxes==ptIdx,:),'butter300T15','butter300TPwidth','nexP50',groupName)
		plotFeatures2D(cellFeatureTable(ptGroupIdxes==ptIdx,:),'spikeT15WidthButter300Order2','spikeTPWidthButter300Order2',groupName,isInterneuron(ptGroupIdxes==ptIdx))
		%plotAllCombinations(cellFeatureTable(ptGroupIdxes==ptIdx,:),[],groupName)
	end

	end

	%close all
end

if(plotGroups)

[ptSzGroupIdxes,ptNames,szNames]=findgroups(cellFeatureTable.PtID,cellFeatureTable.SzNum);

%%
for ptSzIdx=1:length(unique(ptSzGroupIdxes))

	groupName=sprintf('Patient-%s-Seizure-%d',ptNames{ptSzIdx},szNames(ptSzIdx))
	
	%plotFeatures3D(cellFeatureTable(ptSzGroupIdxes==ptSzIdx,:),'butter300T15','butter300TPwidth','nexP50',groupName)
	%plotFeatures2D(cellFeatureTable(ptSzGroupIdxes==ptSzIdx,:),'butter300T15','butter300TPwidth','nexP50',groupName)
		plotFeatures2D(cellFeatureTable(ptSzGroupIdxes==ptSzIdx,:),'spikeT15WidthButter300Order2','spikeTPWidthButter300Order2',groupName,isInterneuron(ptSzGroupIdxes==ptSzIdx))

	%plotAllCombinations(cellFeatureTable(ptSzGroupIdxes==ptSzIdx,:),[],groupName)
	%close all
end
end
%close all

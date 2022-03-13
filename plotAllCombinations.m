function [] = plotAllCombinations(featureTable,featureIdxes,groupName)
	if(isempty(featureIdxes))
		featureIdxes=3:width(featureTable);
	end

	alreadyDone=zeros(length(featureIdxes));	
	featureNames=featureTable.Properties.VariableNames;
	for f1I=1:length(featureIdxes)
		for f2I=1:length(featureIdxes)
			if(f1I~=f2I && ~alreadyDone(f1I,f2I) && ~alreadyDone(f2I,f1I))
				f1=featureIdxes(f1I);
				f2=featureIdxes(f2I);
				close all
				figure
				for cellNum=1:height(featureTable)
					plot(featureTable{cellNum,f1},featureTable{cellNum,f2},'ko','MarkerSize',1)
					hold on
				end
				featureName1=featureNames{f1};
				featureName2=featureNames{f2};
				xlabel(featureName1)
				ylabel(featureName2)
				alreadyDone(f1I,f2I)=1;
				saveName=sprintf('%s-%svs%s',groupName,featureName1,featureName2);
				title(saveName)
				saveas(gcf,[saveName '.tif'])	
			
			end
		end
	end
	
	
	

function [] = plotAllFeatureCombs()
	disp('Started!')
	numFeatures=25;
	alreadyDone=zeros(numFeatures);
	count=1;
	
	filterType='Butter300-3000';

	inputCellArray=cell(1,3);
	startIdx=4;
	for f1=startIdx:numFeatures
		for f2=startIdx:numFeatures
			if(f1~=f2 && alreadyDone(f1,f2)<1 && alreadyDone(f2,f1)<1)
				%tic
				%plot2DFeatureSpaceAllCells(f1,f2);
				%toc
				inputCellArray{count,1}=f1;
				inputCellArray{count,2}=f2;
				inputCellArray{count,3}=filterType;

				alreadyDone(f1,f2)=1;

				count=count+1
			end
		end
	end
inputCellArray
	save('allCellsAllFeatureCombsPar.mat','inputCellArray')

function [U,sigma,V] = svdAccessRoutine(xaxisValues,A,columnDescrp,xUnitStr,iterationUnitStr,yUnitStr,columnCategories,xValueNames,manualColors,nameToClrDict)
%freqBandName='Alpha';
%Description: given matrix mat file and descriptive labels, saves first 20 singular value modes describing data
plot3d=1;

if(~exist('xaxisValues') || sum(isnan(xaxisValues))>0)
	xaxisValues=1:size(A,1);
end

iterationIDs=1:size(A,2);

if(~exist('columnCategories','var'))
	columnCategories={};
	for i=1:size(A,2)
		columnCategories{end+1}='default';
	end
end

if(~iscell(columnCategories))
	columnCategories=cellstr(num2str(columnCategories(:)));
end

rezeroPCspace=0;

recompute=1;
plotSummary=1;
zeroMeanCol=1; %makes PCA



%categoryNames=unique(columnCategories);
if(exist('nameToClrDict','var'))
    categoryNames=nameToClrDict.Keys;
else
    categoryNames={'default'};
end

numCategories=length(categoryNames);
if(~exist('manualColors') || isnan(manualColors(1)))
	categoryColors=copper(numCategories);
else
	categoryColors=manualColors;
end

if(zeroMeanCol)
	%subtract mean column from each column
	A=A-repmat(mean(A,2),1,size(A,2));
end

if(~exist(sprintf('%s-SVD.mat',columnDescrp)) || recompute == 1)
    disp(sprintf('getting SVD of  %s matrix....',columnDescrp))
    
    maxNumCols=min(1000000,size(A,2));
    [U,S,V]=svd(A(:,1:maxNumCols),'econ');
    %[U,S,V]=svd(A(:,1:end),'econ');
    %numVectorsSave=20;
    numVectorsSave=max(5,size(A,2));
    %{
    U=U(:,1:numVectorsSave);
    S=S(:,1:numVectorsSave);
    V=V(:,1:numVectorsSave);
    %}
	%V=V(1:numVectorsSave,:);
    save(sprintf('%s-SVD.mat',columnDescrp),'U','S','V');

else
    svdData=load(sprintf('%s-SVD.mat',columnDescrp));
    U=svdData.U;
    S=svdData.S;
    V=svdData.V;
end

sigma=diag(S);

if(plotSummary)
	figure
	numComponents=min(size(U,1),4);
	subplot(numComponents+1,1,1)


	%plot(sigma(1:numVectorsSave),'r*')
	plot(sigma,'r*')
	xlabel('Component number')
	ylabel('Singular value')
	title(sprintf('%s: Strength of components',columnDescrp))


	%timeAxis=(1:size(U,1))/decFs;
	for i=1:numComponents
	   subplot(numComponents+1,1,i+1)

	   plot(xaxisValues,U(:,i))
	   %plot(timeAxis,log(U(:,i)))
	   %xlim([1 length(U(:,i))])
	   xlabel(xUnitStr)
	   %xlabel('Freq (Hz)')
	   ylabel(yUnitStr)
	   %ylabel('ln(Power)')
	   title(sprintf('Scalable component %d',i))

	end

    figure
    for i=1:numComponents
	   subplot(numComponents+1,1,i+1)

	   plot(iterationIDs,V(:,i))
	   %plot(timeAxis,log(U(:,i)))
	   %xlim([1 length(U(:,i))])
	   xlabel(iterationUnitStr)
	   %xlabel('Freq (Hz)')
	   ylabel(yUnitStr)
	   %ylabel('ln(Power)')
	   title(sprintf('Component %d strength across %s',i,iterationUnitStr))

	end

    

	%displayCurrentFigure(sprintf('%s_componentForms.tif',columnDescrp))
	%displayCurrentFigure(sprintf('%s_componentForms.fig',columnDescrp))
    %{
	figure
	for i=1:numComponents
		subplot(numComponents,1,i)
		if(iscell(xValueNames) || ~isnan(xValueNames))
			pcolorRowWithText(xValueNames,U(:,i))	
		end
		[varExplainedFrac]=getVarExplained(sigma,i)
		title(sprintf('PC%d (scaling this pattern accounts for %.2f percent data)',i,varExplainedFrac*100))
	end
    %}
	figure; 

	if(rezeroPCspace)
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%rezero PC space
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		newPCspaceZero=[0.03802 -0.0827 -0.2436];
		V(:,1:3)=V(:,1:3)-newPCspaceZero;

		[phi,theta,r]=cart2sph(V(:,1),V(:,2),V(:,3));
		figure
		for colNum=1:size(V,1)
			currCategoryNum=strcmp(unique(columnCategories),columnCategories{colNum});
                        currColor=categoryColors(find(currCategoryNum),:);
			plot(phi(colNum),theta(colNum),'.','Color',currColor,'MarkerSize',10)
			%plot(phi(colNum),r(colNum),'.','Color',currColor,'MarkerSize',25)
			%plot(phi(colNum),,'.','Color',currColor,'MarkerSize',25)
                        hold on
		end
		xlabel('Theta')		
		ylabel('Phi')		
		fds
	end

	if(plot3d)
		%plot(V(:,1),V(:,2),'ko')
		usedClrIdx=zeros(size(unique(columnCategories)));
		figure
		for colNum=1:size(V,1)
			%currCategoryNum=strcmp(unique(columnCategories),columnCategories{colNum});
			%currColor=categoryColors(find(currCategoryNum),:);
            if(exist('nameToClrDict','var'))
                clrIdx=nameToClrDict(columnCategories{colNum});
            else
                clrIdx=1;
            end
			currColor=categoryColors(clrIdx,:);
			usedClrIdx(clrIdx)=1;

			plot3(V(colNum,1),V(colNum,2),V(colNum,3),'.','Color',currColor,'MarkerSize',25*2)
			hold on
		end

		%xlim([prctile(V(:,1),0.5) prctile(V(:,1),99.5)])
		%ylim([prctile(V(:,2),0.5) prctile(V(:,2),99.5)])
		%xlim([prctile(V(:,1),0.5) prctile(V(:,1),99.5)])
		%ylim([prctile(V(:,2),0.5) prctile(V(:,2),99.5)])
		xlim([-0.15 0.25])
		ylim([-0.45 0.25])
		zlim([-0.25 0.35])
		currXlim=xlim;
		currYlim=ylim;
		currZlim=zlim;
		xlabel('PC1')
		ylabel('PC2')
		zlabel('PC3')
		hold on
		plot3(currXlim,[0 0],[0 0],'r--')
		plot3([0 0],currYlim,[0 0],'r--')
		plot3([0 0],[0 0],currZlim,'r--')

		if(rezeroPCspace)
			title(sprintf('%s: Rezeroed scale of components per instance',columnDescrp))
		else
			title(sprintf('%s: Scale of components per instance',columnDescrp))

		end
		%addCustomLegendTibin(categoryColors{logical(usedClrIdx)},categoryNames{logical(usedClrIdx)});
		addCustomLegendTibin(categoryColors(logical(usedClrIdx),:),categoryNames(logical(usedClrIdx)));
		setFigFontTo(28)
		%plot 3 pairwise projections
		PCnum1=1;
		PCnum2=2;
		plotPairwisePCspace
		
		PCnum1=1;
		PCnum2=3;
		plotPairwisePCspace
		
		PCnum1=2;
		PCnum2=3;
		plotPairwisePCspace
	else
		%plot(V(:,1),V(:,2),'ko')
		usedClrIdx=zeros(size(unique(columnCategories)));
                for colNum=1:size(V,1)
                        %currCategoryNum=strcmp(unique(columnCategories),columnCategories{colNum});
                        %currColor=categoryColors(find(currCategoryNum),:);
                        clrIdx=nameToClrDict(columnCategories{colNum});
                        currColor=categoryColors(clrIdx,:);
                        usedClrIdx(clrIdx)=1;
			plot(V(colNum,1),V(colNum,2),'o','Color',currColor,'MarkerSize',5)
                        hold on
                end

                %xlim([prctile(V(:,1),0.5) prctile(V(:,1),99.5)])
                %ylim([prctile(V(:,2),0.5) prctile(V(:,2),99.5)])
                %xlim([prctile(V(:,1),0.5) prctile(V(:,1),99.5)])
                %ylim([prctile(V(:,2),0.5) prctile(V(:,2),99.5)])
                currXlim=xlim;
                currYlim=ylim;
                xlabel('PC1')
                ylabel('PC2')
                hold on
                plot(currXlim,[0 0],'r--')
                plot([0 0],currYlim,'r--')
	
		title(sprintf('%s: Scale of components per instance',columnDescrp))
		addCustomLegendTibin(categoryColors{logical(usedClrIdx)},categoryNames{logical(usedClrIdx)});;

	end	
	%displayCurrentFigure(sprintf('%s_componentStrengths.tif',columnDescrp))
	%displayCurrentFigure(sprintf('%s_componentStrengths.fig',columnDescrp))
end

function [W,H] =nnmfAccessRoutine(xaxisValues,A,columnDescrp,xUnitStr,iterationUnitStr,yUnitStr,columnCategories,xValueNames,manualColors,nameToClrDict)
%freqBandName='Alpha';
%Description: given matrix mat file and descriptive labels, saves first 20 singular value modes describing data
plot3d=1;

%numAssemblies=10;
%numAssemblies=15;
numAssemblies=17;
%numAssemblies=20;
%numAssemblies=25;

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
zeroMeanCol=0; %all values must be non-negative



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

if(~exist(sprintf('%s-NNMF.mat',columnDescrp)) || recompute == 1)
    disp(sprintf('getting NNMF of  %s matrix....',columnDescrp))
    
    maxNumCols=min(1000000,size(A,2));
    [W,H]=nnmf(A(:,1:maxNumCols),numAssemblies);
    save(sprintf('%s-NNMF.mat',columnDescrp),'W','H');
else
    nnmfData=load(sprintf('%s-NNMF.mat',columnDescrp));
    W=nnmfData.W;
    H=nnmfData.H;
end

for i=1:size(H,1)
	[cmX,cmVal]=getCenterOfMass(H(i,:));
	hCOMs(i)=cmX;
end

[~,hComSort]=sort(hCOMs)

W=W(:,hComSort);
H=H(hComSort,:);


if(plotSummary)
	figure
	%numComponents=min(size(W,1),4);
	numComponents=min(size(W,2),numAssemblies);
	subplot(numComponents+1,1,1)

	%timeAxis=(1:size(U,1))/decFs;
	for i=1:numComponents
	   subplot(numComponents,1,i)

	   plot(xaxisValues,W(:,i))
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
	   subplot(numComponents,1,i)

	   plot(iterationIDs,H(i,:))
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

	if(plot3d)
		%plot(V(:,1),V(:,2),'ko')
		usedClrIdx=zeros(size(unique(columnCategories)));
		figure
		for colNum=1:size(H,1)
			%currCategoryNum=strcmp(unique(columnCategories),columnCategories{colNum});
			%currColor=categoryColors(find(currCategoryNum),:);
		    if(exist('nameToClrDict','var'))
			clrIdx=nameToClrDict(columnCategories{colNum});
		    else
			clrIdx=1;
		    end
			currColor=categoryColors(clrIdx,:);
			usedClrIdx(clrIdx)=1;

			plot3(H(1,colNum),H(2,colNum),H(3,colNum),'.','Color',currColor,'MarkerSize',25*2)
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
		%{
		PCnum1=1;
		PCnum2=2;
		plotPairwisePCspace
		
		PCnum1=1;
		PCnum2=3;
		plotPairwisePCspace
		
		PCnum1=2;
		PCnum2=3;
		plotPairwisePCspace
		%}
	else
		%plot(V(:,1),V(:,2),'ko')
		usedClrIdx=zeros(size(unique(columnCategories)));
                for colNum=1:size(H,1)
                        %currCategoryNum=strcmp(unique(columnCategories),columnCategories{colNum});
                        %currColor=categoryColors(find(currCategoryNum),:);
			    if(exist('nameToClrDict','var'))
				clrIdx=nameToClrDict(columnCategories{colNum});
			    else
				clrIdx=1;
			    end
                        currColor=categoryColors(clrIdx,:);
                        usedClrIdx(clrIdx)=1;
			plot(H(1,colNum),H(2,colNum),'o','Color',currColor,'MarkerSize',5)
                        hold on
                end

                %xlim([prctile(V(:,1),0.5) prctile(V(:,1),99.5)])
                %ylim([prctile(V(:,2),0.5) prctile(V(:,2),99.5)])
                %xlim([prctile(V(:,1),0.5) prctile(V(:,1),99.5)])
                %ylim([prctile(V(:,2),0.5) prctile(V(:,2),99.5)])
                currXlim=xlim;
                currYlim=ylim;
                xlabel('W1')
                ylabel('W2')
                hold on
                plot(currXlim,[0 0],'r--')
                plot([0 0],currYlim,'r--')
	
		title(sprintf('%s: Scale of components per instance',columnDescrp))
		addCustomLegendTibin(categoryColors{logical(usedClrIdx)},categoryNames{logical(usedClrIdx)});;

	end	
	%displayCurrentFigure(sprintf('%s_componentStrengths.tif',columnDescrp))
	%displayCurrentFigure(sprintf('%s_componentStrengths.fig',columnDescrp))
end

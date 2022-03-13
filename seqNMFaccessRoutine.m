function [W,H] =seqNMFaccessRoutine(xaxisValues,A,columnDescrp,xUnitStr,iterationUnitStr,yUnitStr,figComps,figIters,columnCategories,xValueNames,manualColors,nameToClrDict)
%freqBandName='Alpha';
%Description: given matrix mat file and descriptive labels, saves first 20 singular value modes describing data
plot3d=0;

%numAssemblies=10;
%numAssemblies=15;
%numAssemblies=17;
%numAssemblies=20;
%numAssemblies=25;
numAssemblies=50;

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

if(~exist(sprintf('%s-seqNMF.mat',columnDescrp)) || recompute == 1)
    disp(sprintf('getting seqNMF of  %s matrix....',columnDescrp))
    
    maxNumCols=min(1000000,size(A,2));
    %[W,H]=nnmf(A(:,1:maxNumCols),numAssemblies);
    %K=number of factors
    %L=convolutional width
    lambda =.005; %step size
    %computing seqNMF
    [W,H]=seqNMF(A(:,1:maxNumCols),'K',numAssemblies,'L',1,'lambda',lambda);
    save(sprintf('%s-seqNMF.mat',columnDescrp),'W','H');
else
    nnmfData=load(sprintf('%s-seqNMF.mat',columnDescrp));
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
	figure(figComps)
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

    figure(figIters)
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
end

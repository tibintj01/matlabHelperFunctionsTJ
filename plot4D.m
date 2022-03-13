function [] = plot4D(x,y,z,c,numBins,cEdges)
%UNTITLED Summary of this function goes here
% c is a 3D matrix with domain given by x,y,z

%maxMarkerWidth=70;
maxMarkerWidth=50;
maxMarkerWidth=60*(10/numBins);
lineWidth=5;
lineWidth=2;

xEdges=linspace(min(x),max(x),numBins+1);
yEdges=linspace(min(y),max(y),numBins+1);
zEdges=linspace(min(z),max(z),numBins+1);

if(~exist('cEdges','var'))
    cEdges=linspace(min(c(:)),max(c(:)),numBins+1);
else
    cEdges=linspace(cEdges(1),cEdges(end),numBins+1);
end

c(c>cEdges(end))=cEdges(end);
c(c<cEdges(1))=cEdges(1);


xBins=edgesToBins(xEdges);
yBins=edgesToBins(yEdges);
zBins=edgesToBins(zEdges);
cBins=edgesToBins(cEdges);

xInBins=discretize(x,xEdges);
yInBins=discretize(y,yEdges);
zInBins=discretize(z,zEdges);
cInBins=discretize(c,cEdges);

colors=jet(length(cBins));
%colors=parula(length(cBins));
figure


for xBin=1:length(xBins)
    disp(xBin/length(xBins))
	for yBin=1:length(yBins)
		for zBin=1:length(zBins)
			curr3DptColor=cInBins(xBin,yBin,zBin);
            if(isnan(curr3DptColor))
                continue
            end
			plot(xInBins(xBin),yInBins(yBin),'o','MarkerSize',zInBins(zBin)/max(zInBins)*maxMarkerWidth,'Color',colors(curr3DptColor,:),'LineWidth',lineWidth)
			hold on
		end
	end
end
xlim([xEdges(1)-0.5 xEdges(end)+0.5])
ylim([yEdges(1)-0.5 yEdges(end)+0.5])

daspect([1 1 1])



N=1000;

%numNull=100000;
numNull=1000;

nullNormDIvals=NaN(numNull,1);

maxDI=getDI(N:-1:1);

for i=1:numNull
	i
	randDI=getDI(randperm(N));
	
	nullNormDIvals(i)=randDI/maxDI;	
end

speed1=0.25;
speed2=1.25;

speed1Vals=[0.5255
    0.5355
    0.5198
    0.5223
    0.5236
    0.5296
    0.5249
    0.5225
    0.5254
    0.5232];

speed2Vals=[  0.5309
    0.5223
    0.5198
    0.5411
    0.5454
    0.5318
    0.5289
    0.5215
    0.5276
    0.5304];

%edges=0:0.01:1;
edges=0:0.01:1;
binCenters=edgesToBins(edges);

meanSpeed1=mean(speed1Vals);

meanSpeed2=mean(speed2Vals);


[nullDIdist]=histcounts(nullNormDIvals,edges);

[speed1DIdist]=histcounts(speed1Vals,edges);
[speed2DIdist]=histcounts(speed2Vals,edges);

figure
plot(binCenters,nullDIdist/sum(nullDIdist),'k','LineWidth',3)
hold on

plot(binCenters,speed1DIdist/sum(speed1DIdist),'b','LineWidth',3)
plot(binCenters,speed2DIdist/sum(speed2DIdist),'r','LineWidth',3)

yMax=0.8;
plot([meanSpeed1 meanSpeed1],[0 yMax],'b--','LineWidth',1)
plot([meanSpeed2 meanSpeed2],[0 yMax],'r--','LineWidth',1)
%xlim([0 1])
xlim([0.44 0.58])
ylim([0 yMax])

xlabel('normalized DI value')
ylabel('Probability density')

legend(sprintf('random permutation DI (n=%d)',numNull), sprintf('speed=%.2f network output',speed1), sprintf('speed=%.2f network output',speed2),'Location','Best')

setFigFontTo(18)

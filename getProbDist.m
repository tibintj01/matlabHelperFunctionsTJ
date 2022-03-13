function [p] = getProbDist(vals,edges,showPlots,isCirc)

smoothWindBins=10;

[N,c]=histcounts(vals,edges);
bins=edgesToBins(edges);

p=N/sum(N(:));

numBins=length(p);

if(isCirc)
    pRepeated=[p(:); p(:); p(:)];
    pRepeated=smooth(pRepeated,smoothWindBins);
    
    startBin=numBins+1;
    endBin=startBin+(numBins-1);
    
    p=pRepeated(startBin:endBin);
else
    p=smooth(p,smoothWindBins);
end

p=p/sum(p(:));
if(showPlots)
    %figure

    yyaxis left
    plot(bins,N,'LineWidth',1)

    yyaxis right
    plot(bins,p,'LineWidth',5)
end





function [pL] = plotLineHist(vals,edges,colorStr,plotZeroLine)
%plots histogram as line plot rather than bar plot

if(~exist('colorStr','var'))
    colorStr='k-';
end
    
    if(~exist('plotZeroLine','var'))
        plotZeroLine=0;
    end
    
    [N]=histcounts(vals,edges);
    Nsmoothed=smooth(N)/sum(smooth(N));
    hold on
    pL=plot(edgesToBins(edges),Nsmoothed,colorStr,'LineWidth',5)
    
    if(plotZeroLine)
        hold on
        plot([0 0], ylim,'k--','LineWidth',5)
    end
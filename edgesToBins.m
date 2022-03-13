function [bins] = edgesToBins(edges)
	binWidth=median(diff(edges));
	
	bins=edges+binWidth/2;
	bins=bins(1:(end-1));

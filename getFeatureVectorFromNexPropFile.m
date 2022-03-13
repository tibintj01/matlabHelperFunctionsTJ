function [featureVector] = getFeatureVectorFromFile(filePath)
	ptID=getPtID(filePath);
	szID=getSzID(filePath);
	butterFilePath=strrep(filePath,'nex','butter300');
	rawFilePath=strrep(filePath,'nex','raw');

	cellPropsNex=load(filePath);

	%get good spike indices from nex-filter setting properties
	goodSpikeIdxes=cellPropsNex.goodSpikes;

    %if(~exist(butterFilePath,'file') || ~exist(rawFilePath,'file'))
    %    featureVector={'NaN',NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN};
    %    return
    %end
    
	cellPropsButter=load(butterFilePath);
	cellPropsRaw=load(rawFilePath);

	
	nexTPwidth=nanmean(cellPropsNex.spikeOverallMaxTPWidth(goodSpikeIdxes));
	butter300TPwidth=nanmean(cellPropsButter.spikeOverallMaxTPWidth(goodSpikeIdxes));
	rawTPwidth=nanmean(cellPropsRaw.spikeOverallMaxTPWidth(goodSpikeIdxes));

	nexP50=nanmean(cellPropsNex.spikeP50Width(goodSpikeIdxes));
	nexP25=nanmean(cellPropsNex.spikeP25Width(goodSpikeIdxes));
	nexT50=nanmean(cellPropsNex.spikeT50Width(goodSpikeIdxes));
	

	rawT15=nanmean(cellPropsRaw.spikeT15Width(goodSpikeIdxes));
	butter300T15=nanmean(cellPropsButter.spikeT15Width(goodSpikeIdxes));
	nexT15=nanmean(cellPropsNex.spikeT15Width(goodSpikeIdxes));

	butter300P25=nanmean(cellPropsButter.spikeP25Width(goodSpikeIdxes));
	butter300T50=nanmean(cellPropsButter.spikeT50Width(goodSpikeIdxes));

	featureVector={ptID,szID,nexTPwidth,rawTPwidth,butter300TPwidth,nexP50,nexP25,nexT50,rawT15,butter300T15,nexT15,butter300P25,butter300T50};
end

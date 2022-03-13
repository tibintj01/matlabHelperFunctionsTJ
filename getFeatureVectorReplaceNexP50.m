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


	nexP50=nanmean(cellPropsNex.spikeP50Width(goodSpikeIdxes));
	butter300TPwidth=nanmean(cellPropsButter.spikeOverallMaxTPWidth(goodSpikeIdxes));
	butter300T15=nanmean(cellPropsButter.spikeT15Width(goodSpikeIdxes));
	butter300T30=nanmean(cellPropsButter.spikeT30Width(goodSpikeIdxes));	
	butter300T75=nanmean(cellPropsButter.spikeT75Width(goodSpikeIdxes));
	butter300TPwidthNotOverall=nanmean(cellPropsButter.spikeTPWidth(goodSpikeIdxes));

	butter300P25=nanmean(cellPropsButter.spikeP25Width(goodSpikeIdxes));
	butter300T50=nanmean(cellPropsButter.spikeT50Width(goodSpikeIdxes));
	butter300P75=nanmean(cellPropsButter.spikeT75Width(goodSpikeIdxes));	

	featureVector={ptID,szID,nexP50,butter300TPwidth,butter300T15,butter300T30,butter300T75,butter300TPwidthNotOverall,butter300P25,butter300T50,butter300P75};

end

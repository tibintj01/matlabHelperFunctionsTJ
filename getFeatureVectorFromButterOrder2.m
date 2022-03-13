function [featureVals,featureNames] = getFeatureVectorFromFile(filePath,useAllSpikes)
	try
		ptID=getPtID(filePath);
		szID=getSzID(filePath);
	catch
		disp('Sz session not detected...')
		ptID=0;
		szID=0;
	end

	cellPropsButter=load(filePath);

	%get good spike indices from nex-filter setting properties
	goodSpikeIdxes=cellPropsButter.goodSpikes;

    %if(~exist(butterFilePath,'file') || ~exist(rawFilePath,'file'))
    %    featureVector={'NaN',NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN};
    %    return
    %end
   	fieldNames=fieldnames(cellPropsButter);

	featureVals=cell(33,1);
	featureVals{1}=ptID;
	featureVals{2}=szID;
	featureNames={'PtID','SzNum'};
	for i=5:34
		featureArray=cellPropsButter.(fieldNames{i});
		if(exist('useAllSpikes','var') && useAllSpikes)
			featureVals{i-2}=nanmean(featureArray(:));
		else
			featureVals{i-2}=nanmean(featureArray(goodSpikeIdxes));
		end
		featureNames{end+1}=[fieldNames{i} 'Butter300Order2'];
	end
	

	troughToPeakRatioButter=cellPropsButter.spikeMinAmp(goodSpikeIdxes)./cellPropsButter.spikeMaxAmp(goodSpikeIdxes);

	featureNames=[featureNames 'troughToPeakRatioButter300Order2'];
	featureVals{33}=nanmean(troughToPeakRatioButter);
	featureVals=featureVals';
end

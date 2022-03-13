function [var1Bins,var2Bins,spaceTimeRateMeansSmooth,dataCountPerBin] = makeHeatMapAvgOf3rdVar(var1,var2,var3,numBins)

    var1=var1(:);
    var2=var2(:);
    var3=var3(:);
    
    nanIDs=isnan(var1(:)) | isnan(var2(:)) | isnan(var3(:));
    var1(nanIDs)=[];
    var2(nanIDs)=[];
    var3(nanIDs)=[];

    [var1Disc,var1Edges]=discretize(var1(:),numBins);
    [var2Disc,var2Edges]=discretize(var2(:),numBins);
    
    var1Bins=edgesToBins(var1Edges);
    var2Bins=edgesToBins(var2Edges);


    data=[var1Disc(:)'; var2Disc(:)'; var3(:)'];
    
    size(data)
    [zs,avgs] = threeDPtsToSurface(data);

    dataCountPerBin=zeros(numBins,numBins);    
    spaceTimeRateMeans=NaN(size(zs));
    for i=1:numBins
        for j=1:numBins
           if(length(zs{i,j})>0)
               spaceTimeRateMeans(i,j)=nanmean(zs{i,j});
		currBinData=zs{i,j};
		currBinData=currBinData(~isnan(currBinData));
		dataCountPerBin(i,j)=dataCountPerBin(i,j)+length(currBinData);
           end
        end
    end
    
    filtWidth = 3;
filtSigma = 1.5;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
spaceTimeRateMeansSmooth = nanconv(spaceTimeRateMeans,imageFilter, 'nanout','edge');
    
  

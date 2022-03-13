function [xCenters,yCenters,spaceTimeRateMeansSmooth,dataCountPerBin] =makeHeatMapOf3rdVarUniv(input)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %input struct contains data in 3 cols, x and y radii, x and y moving center edges
    %and desired statistic for locally collected values. Required field names:
    %	input.data3cols
    %   input.xRad
    %	input.yRad
    %	input.xEdges
    %	input.yEdges
    %	input.desiredStat
    %Desired stat options:
    %	'mean'
    %   'SEM'
    %   'circMean'
    %   'circPeakPhase'
    %   'circR'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data3cols=input.data3cols;
    xRad=input.xRad;
    yRad=input.yRad;
    xEdges=input.xEdges;
    yEdges=input.yEdges;
    desiredStat=input.desiredStat;

    var1=data3cols(:,1);
    var2=data3cols(:,2);
    var3=data3cols(:,3);
    
    nanIDs=isnan(var1(:)) | isnan(var2(:)) | isnan(var3(:));
    var1(nanIDs)=[];
    var2(nanIDs)=[];
    var3(nanIDs)=[];

    data3cols=[var1(:), var2(:), var3(:)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %old approach: discretize x and y, accumulate corresponding z's
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    [var1Disc,var1Edges]=discretize(var1(:),numBins);
    [var2Disc,var2Edges]=discretize(var2(:),numBins);
    var1Bins=edgesToBins(var1Edges);
    var2Bins=edgesToBins(var2Edges);
    data=[var1Disc(:)'; var2Disc(:)'; var3(:)'];
    [zs,avgs] = threeDPtsToSurface(data);
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %new approach: moving window local  x and y, accumulate corresponding z's
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %zs=[];
      
      
    xCenters=edgesToBins(xEdges);
    yCenters=edgesToBins(yEdges);

    numXcenters=length(xCenters);
    numYcenters=length(yCenters);

    zs=cell(length(xCenters),length(yCenters));
    
    for ci=1:length(xCenters)
        disp(sprintf('calculating average for row %d out of %d', ci,length(xCenters)))
        tic
        for cj=1:length(yCenters)
            currXcenter=xCenters(ci);
            currYcenter=yCenters(cj);

            currCenter=[currXcenter; currYcenter];
            [currWindowData]=getDataWithinRect(data3cols,currCenter,xRad,yRad);
            if(~isnan(currWindowData))
                zs{ci,cj}=currWindowData(:,3);
            else
                zs{ci,cj}=NaN;
            end
        end
        toc
    end 



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get average or desired statistic on each locally collected value set
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataCountPerBin=zeros(numXcenters,numYcenters);    
    spaceTimeRateMeans=NaN(size(zs));
    for i=1:numXcenters
        for j=1:numYcenters
            if(i>size(zs,1))
                continue
            end
            
            if(j>size(zs,2))
                continue
            end
            
           %if(length(zs{i,j})>0 && (sum(isnan(zs{i,j})) ~= length(zs{i,j})))
           if(length(zs{i,j})>1 && (sum(isnan(zs{i,j})) ~= length(zs{i,j}))) %ignore single points in window
                currBinData=zs{i,j};
                
		if(strcmp(desiredStat,'mean'))
		    spaceTimeRateMeans(i,j)=nanmean(currBinData);
		elseif(strcmp(desiredStat,'fractionSEM'))
                    spaceTimeRateMeans(i,j)=getSEMacrossRows(currBinData(:))/nanmean(currBinData(:));
		elseif(strcmp(desiredStat,'SEM'))
                    spaceTimeRateMeans(i,j)=getSEMacrossRows(currBinData(:));
		elseif(strcmp(desiredStat,'circMean'))
               	    spaceTimeRateMeans(i,j)=circMeanDeg(currBinData(:));
		elseif(strcmp(desiredStat,'circR'))
                    spaceTimeRateMeans(i,j)=circ_r(ang2rad(currBinData(:)));
        elseif(strcmp(desiredStat,'circPeakPhase'))
                [~,currBinPeakPhase] = getCircKernelDistr(currBinData(:));
                spaceTimeRateMeans(i,j)=currBinPeakPhase;
        end

		currBinData=currBinData(~isnan(currBinData));
		dataCountPerBin(i,j)=dataCountPerBin(i,j)+length(currBinData);
           end
        end
    end
    
    filtWidth = 3;
filtSigma = 1.5;
%filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
spaceTimeRateMeansSmooth = nanconv(spaceTimeRateMeans,imageFilter, 'nanout','edge');
    %spaceTimeRateMeansSmooth=spaceTimeRateMeans;
  

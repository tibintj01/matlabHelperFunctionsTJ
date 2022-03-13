close all; clear all; clc
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData';
saveDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages';

filePaths=getFilePathsRegex(dataDir,'*Info.mat');

for fi=1:length(filePaths)
    currFilePath=filePaths{fi};
    
    currUnitInfo=load(currFilePath);
    currUnitInfo=currUnitInfo.unitInfo
    
    mazeTaskTimeBounds=currUnitInfo.mazeTaskTimeBounds;
    
    firingRateOverTime=currUnitInfo.firingRateOverTime;
    rateTimeBinCenters=currUnitInfo.rateTimeBinCenters;
    
    positionStruct=currUnitInfo.positionPerTime;
    
    currSessionName=currUnitInfo.sessionName;
    
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute position along track from 2D after filling in NaN (LED out of view)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    positionPerTimeStep=positionStruct.OneDLocation; %in METERS
    positionPerTimeStepFilled = fillmissing(positionPerTimeStep,'linear')
    positionPerTimeStepFilled(positionPerTimeStepFilled<0)=NaN;
    %}
    
    positionTimeAxis=positionStruct.TimeStamps;
    approxPosTimeStep=median(diff(positionTimeAxis));
    
    %linear tracks are horizontal - just take 1st col of 2D location 
    if(contains(positionStruct.MazeType,'Linear Maze'))
        positionPerTimeStep=positionStruct.TwoDLocation(:,1);
        %continue
    end
    
  
    %circular track - project onto circle of center (0.5,0.5) and radius 0.5
    if(strcmp('Circular Maze',positionStruct.MazeType)) %slightly less than pi meters circumference
        
        delXPerTimeStep=[diff(positionStruct.TwoDLocation(:,1)); NaN];
        delYPerTimeStep=[diff(positionStruct.TwoDLocation(:,2)); NaN];
        delTPerTimeStep=[diff(positionStruct.TimeStamps(:)); NaN];
        
        distPerTimeStep=sqrt(delXPerTimeStep.^2 + delYPerTimeStep.^2);
        estRawSpeedPerTimeStep=distPerTimeStep./delTPerTimeStep;
        
        %lingeringIDs=estRawSpeedPerTimeStep<0.2; %m/s
         lingeringIDs=estRawSpeedPerTimeStep<0.3; %m/s
          if(strcmp(currSessionName,'Gatsby_08282013'))
              lingeringIDs=estRawSpeedPerTimeStep<0.5; %m/s
          end
        binFact=2;
        [N,c]=hist3(positionStruct.TwoDLocation(~lingeringIDs,:),binFact*[100 100]);
        rowBins=c{1};
        colBins=c{2};
        rowBins=rowBins(:);
        colBins=colBins(:);
        histImgSize=size(N);
        
        N=imclose(N,ones(3));
        
        countThresh=10/binFact;
        N(N<countThresh)=0;
        N(N>=countThresh)=1;
        
        N=medfilt2(N);
        N(N>0)=1;
        
        N = bwareaopen(N,100);
        
        [denseBinRows,denseBinCols]=ind2sub(histImgSize,find(N));
        denseBinCoords=[rowBins(denseBinRows) colBins(denseBinCols)];
        
        %figure; imagesc(N); colorbar
        
        
        %[x0 y0 R]=circle_fit(positionStruct.TwoDLocation(:,1),positionStruct.TwoDLocation(:,2));
        %[x0 y0 R]=circle_fit(positionStruct.TwoDLocation(~lingeringIDs,1),positionStruct.TwoDLocation(~lingeringIDs,2));
        [x0 y0 R]=circle_fit(denseBinCoords(:,1),denseBinCoords(:,2));
       

                
             
                
         
     
       
        centered2Dlocation=[positionStruct.TwoDLocation(:,1)-x0 positionStruct.TwoDLocation(:,2)-y0];
        [positionPerTimeStep,~] = cart2pol(centered2Dlocation(:,1), centered2Dlocation(:,2));
	     %minPosition=min(positionStruct.OneDLocation);
	     %maxPosition=max(positionStruct.OneDLocation);
            minPosition=0;
	     maxPosition=2*pi*R;
         
	    positionPerTimeStep=scaledata(positionPerTimeStep,minPosition,maxPosition);
        
        %{
         %if(strcmp(currSessionName,'Cicero_09102014'))
         if(strcmp(currSessionName,'Gatsby_08282013'))
         %if(strcmp(currSessionName,'Achilles_11012013'))
                figure; 
                subplot(1,2,1)
                plot(denseBinCoords(:,1),denseBinCoords(:,2),'k.');
                daspect([1 1 1])
                subplot(1,2,2)
                plot(positionStruct.TwoDLocation(~lingeringIDs,1),positionStruct.TwoDLocation(~lingeringIDs,2),'k.'); hold on; plot(xlim,[0.5 0.5],'k--'); plot(xlim,[-0.5 -0.5],'k--'); hold on; plot([0.5 0.5],ylim,'k--'); plot([-0.5 -0.5],ylim,'k--');plot([0 0],ylim,'k--'); plot(xlim,[0 0],'k--'); daspect([1 1 1])
                viscircles([x0 y0],R)
                
             figure; plot(positionStruct.TimeStamps(:),positionPerTimeStep,'k-')
            disp('')
          end
        %}
        
        %continue
    end
    
    %{
    figure; 
    subplot(2,1,1)
    plot(positionTimeAxis,positionPerTimeStep)
    
    subplot(2,1,2)
    maxDispIdx=100000;
    minDispIdx=10000;
    %{
    for di=minDispIdx:maxDispIdx
        plot(positionStruct.TwoDLocation(minDispIdx:di,1),positionStruct.TwoDLocation(minDispIdx:di,2))
        drawnow
    end
    %}
    
     plot(positionStruct.TwoDLocation(:,1),positionStruct.TwoDLocation(:,2))
    %}

    
    
    polynomialFitOrder=1; %fit lines over acc time width
    speedEstTimeWidth=0.125; %sec %~5 point linear fit estimation of speed
      %speedEstTimeWidth=0.075; %sec %~5 point linear fit estimation of speed
    supportWindLength=ceil(speedEstTimeWidth/approxPosTimeStep);
    speedPerTimeStep=movingslope(positionPerTimeStep,supportWindLength,polynomialFitOrder,approxPosTimeStep);
    
    absSpeedPerTimeStep=abs(speedPerTimeStep); %meters/sec
    
    %speedPerTimeStep=[diff(positionPerTimeStep); NaN]/approxTimeStep;
    
    %{
    figure; histogram(speedPerTimeStep)
    figure;
    yyaxis left
    plot(rateTimeBinCenters,firingRateOverTime)
    yyaxis right
    plot(positionStruct.TimeStamps,positionPerTimeStep,'k-')
    hold on
        %plot(positionStruct.TimeStamps,positionPerTimeStepFilled,'ko')
    plot(positionStruct.TimeStamps,speedPerTimeStep)
    xlim(mazeTaskTimeBounds)
    %}
    %disp('')
    
    
    save(currFilePath,'positionPerTimeStep','speedPerTimeStep','absSpeedPerTimeStep','positionTimeAxis','approxPosTimeStep','mazeTaskTimeBounds','-append')
end

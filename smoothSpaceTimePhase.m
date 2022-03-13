function [smoothedMatrix] = smoothSpaceTimePhase(spaceTimePhaseMatrix,smoothSpaceSize,smoothTimeSize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
setTightSubplots_SpaceTime
    plotSmoothed=0;
    newN=100;
    sigma=10;
    sigma=15;
     
     %sigma=25;
      %sigma=30;
    %sigma=7.5;
    sigma=10;
    sigma=20;
    originalNumSpacePts=size(spaceTimePhaseMatrix,1);
    originalNumTimePts=size(spaceTimePhaseMatrix,2);
    
    spaceBins=1:originalNumSpacePts;
    timeBins=1:originalNumTimePts;
    
    %[X,Y]=meshgrid(spaceBins,timeBins);
    [X,Y]=meshgrid(timeBins,spaceBins);
    
    
    V=spaceTimePhaseMatrix;
    
    newSpaceBins=linspace(1,originalNumSpacePts,newN);
    newTimeBins=linspace(1,originalNumTimePts,newN);
    
    %[Xq,Yq]=meshgrid(newSpaceBins,newTimeBins);
    [Xq,Yq]=meshgrid(newTimeBins,newSpaceBins);
    
    resampledMatrix=interp2(X,Y,V,Xq,Yq);
    
    
    %K=gaussian_filter(sigma,sigma); %var vs n?
    smoothSizeTime=sigma;
            smoothSizeSpace=sigma;
        K=1/(smoothSizeSpace*smoothSizeTime)*(ones(smoothSizeSpace,smoothSizeTime));
    smoothedMatrixHiRes=imgaussfilt(resampledMatrix,sigma,'FilterDomain','spatial');
    %smoothedMatrixHiRes=nanconv(resampledMatrix,K,'edge');
    smoothedMatrix=interp2(Xq,Yq,smoothedMatrixHiRes,X,Y);
    %smoothedMatrix=smoothedMatrix';
    
    
    
    if(plotSmoothed)
    fH=figure;
    subplot(2,1,1)
    omarPcolor(newSpaceBins,newTimeBins,resampledMatrix',fH)
    colormap(jet)
    
      subplot(2,1,2)
    omarPcolor(newSpaceBins,newTimeBins,smoothedMatrix',fH)
    
     colormap(jet)
    end

end


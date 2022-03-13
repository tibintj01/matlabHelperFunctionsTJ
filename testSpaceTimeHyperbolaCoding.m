if(~exist('allStrongFieldDistPhaseTimeZTriplets.mat','file'))
    
dataDirStrong='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong';
dataDirMaybe='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/maybe';
filePaths=getFilePathsRegex(dataDirStrong,'*mat');
filePathsMaybe=getFilePathsRegex(dataDirMaybe,'*mat');
%filePaths=[filePaths filePathsMaybe];
 
close all
allFieldPhases=[];
allFieldDists=[];
allFieldTimesZ=[];
allFieldSpeeds=[];
for fi=1:length(filePaths)
    currFile=filePaths{fi};
    fi
    data=load(currFile);
    
    normPositions=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,1};
    elapsedTimes=data.spaceTimePhaseInfo.spaceTimePhasePerCycle{:,2};
    
    maxDT=max(normPositions.*elapsedTimes);
    %fH=figure;
    %subplot(2,1,1)
    %plot(normPositions*120,maxDT-normPositions.*elapsedTimes,'k.')
    
    %figure
    plotPhaseVsDistanceAndTime
    allFieldPhases=[allFieldPhases;allPhases(:)];
    allFieldDists=[allFieldDists;allDists(:)];
    allFieldTimesZ=[allFieldTimesZ;allTimes(:)];
    allFieldSpeeds=[allFieldSpeeds;allSpeeds(:)];
end

%tripletTable=table(allFieldDists(:),allFieldTimesZ(:),allFieldPhases(:));
tripletTable=table(allFieldDists(:),allFieldSpeeds(:),allFieldPhases(:));
tripletMatrix=table2array(tripletTable);

save('allStrongFieldDistPhaseTimeZTriplets.mat','tripletTable','tripletMatrix')
else
    load('allStrongFieldDistPhaseTimeZTriplets.mat')
end
%spaceEdges=linspace(0,1,31);
spaceEdges=linspace(0,1,21);
timeEdges=linspace(15,42,26);
spaceBinWidth=spaceEdges(2)-spaceEdges(1);
timeBinWidth=timeEdges(2)-timeEdges(1);

spaceBins=edgesToBins(spaceEdges);
timeBins=edgesToBins(timeEdges);
numSpaceBins=length(spaceBins);
numTimeBins=length(timeBins);

spaceTimePhaseMap=NaN(numSpaceBins,numTimeBins);
for si=1:length(spaceBins)
    for ti=1:length(timeBins)
        currSpaceBinStart=spaceBins(si)-spaceBinWidth/2;
        currTimeBinStart=timeBins(ti)-timeBinWidth/2;
          currSpaceBinEnd=spaceBins(si)+spaceBinWidth/2;
        currTimeBinEnd=timeBins(ti)+timeBinWidth/2;
        
        inCurrTimeBinIdxes=tripletMatrix(:,2)>=currTimeBinStart & tripletMatrix(:,2)<currTimeBinEnd;
        inCurrSpaceBinIdxes=tripletMatrix(:,1)>=currSpaceBinStart & tripletMatrix(:,1)<currSpaceBinEnd;
        
        inCurrBinIdxes=inCurrTimeBinIdxes & inCurrSpaceBinIdxes;
        currBinPhases=tripletMatrix(inCurrBinIdxes,3);
        currBinMeanPhase=circMeanDeg(currBinPhases);
        spaceTimePhaseMap(si,ti)=currBinMeanPhase;
    end
end

%spaceTimePhaseMap=imgaussfilt(spaceTimePhaseMap,0.75);
%filtWidth = 7;
%filtSigma = 5;
filtWidth = 3;
filtSigma = 1.5;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
spaceTimePhaseMap = nanconv(spaceTimePhaseMap,imageFilter, 'nanout','edge');

omarPcolor(spaceBins,timeBins,spaceTimePhaseMap')
%contour(spaceBins,timeBins,spaceTimePhaseMap')
colormap(jet)
cb=colorbar
%caxis([0 360])
xlabel('Distance in field (fraction)')
ylabel('Running speed (cm/s)')
ylabel(cb,'Avg theta phase (deg)')

title({'Phase precession distance-speed surface','N=94 strongly precessing place cells'})
setFigFontTo(18)
maxFig
saveas(gcf,'PhasePrecessionDistanceSpeedSurface.tif')

dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong';
imageDir='/Users/tibinjohn/thetaSeq/code/spacetimeBehaviorAdjustFieldBounds';
close all
%filePaths=getFilePathsRegex(imageDir,'*tif');
filePaths=getFilePathsRegex(imageDir,'*png');

for fi=1:length(filePaths)
    currFilePath=filePaths{fi};
    
    fileName=getFileNameFromPath(currFilePath);
   
    
    imshow(currFilePath)
    
    if(fi==1)
     maxFig
    end
    drawnow
    
    prompt='Place field/precession start (cm)?';
    manualFieldStartCmAdjusted=input(prompt);
    
    prompt='Place field/precession end (cm)?';
    manualFieldEndCmAdjusted=input(prompt);
 
     
    %dataFileName=[fileName(1:(end-4)) '_spaceTimePhaseCycleInfo.mat'];
    dataFileName=[fileName(20:(end-4)) '_spaceTimePhaseCycleInfo.mat'];
    dataFilePath=fullfile(dataDir,dataFileName);
    
    save(dataFilePath,'manualFieldStartCmAdjusted','manualFieldEndCmAdjusted','-append')
    data=load(dataFilePath)
end

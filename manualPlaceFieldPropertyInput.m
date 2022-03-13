imageDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong/tifs';
dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/manuallySortedPlaceCellsAbove10Hz/strong';


filePaths=getFilePathsRegex(imageDir,'*tif');

for fi=1:length(filePaths)
    currFilePath=filePaths{fi};
    
    fileName=getFileNameFromPath(currFilePath);
   
    
    imshow(currFilePath)
    
    if(fi==1)
     maxFig
    end
    drawnow
    
    prompt='Place field/precession start (cm)?';
    manualFieldStartCm=input(prompt);
    
    prompt='Place field/precession end (cm)?';
    manualFieldEndCm=input(prompt);
 
    prompt='Phase precession start (degrees)?';
    manualPrecessionStartPhaseDeg=input(prompt);
     
    dataFileName=[fileName(1:(end-4)) '_spaceTimePhaseCycleInfo.mat'];
    dataFilePath=fullfile(dataDir,dataFileName);
    
    save(dataFilePath,'manualFieldStartCm','manualFieldEndCm','manualPrecessionStartPhaseDeg','-append')
    data=load(dataFilePath)
end
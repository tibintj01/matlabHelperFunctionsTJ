imageDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/allPlaceCellsPeak10HzWholeTrack';
dataDir='/Users/tibinjohn/thetaSeq/code/spaceTimeProcessedDir/';

manualPlaceCellFail=fullfile(dataDir,'manuallySortedPlaceCellsAbove10Hz','fail');
manualPlaceCellStrong=fullfile(dataDir,'manuallySortedPlaceCellsAbove10Hz','strong');
manualPlaceCellMaybe=fullfile(dataDir,'manuallySortedPlaceCellsAbove10Hz','maybe');

touchDir(manualPlaceCellFail);
touchDir(manualPlaceCellStrong);
touchDir(manualPlaceCellMaybe);

filePaths=getFilePathsRegex(imageDir,'*tif');

for fi=1:length(filePaths)
    currFilePath=filePaths{fi};
    
    fileName=getFileNameFromPath(currFilePath);
   
    
    imshow(currFilePath)
    if(fi==1)
     maxFig
    end
    drawnow
    
    prompt='Place cell precession quality? f=fail, s=strong, ELSE=maybe';
    q=input(prompt,'s');
 
     
    dataFileName=[fileName(1:(end-4)) '_spaceTimePhaseCycleInfo.mat'];
    dataFilePath=fullfile(dataDir,dataFileName);
    
    if(q=='f')
        movefile(dataFilePath, manualPlaceCellFail)
        movefile(currFilePath, manualPlaceCellFail)
        disp('fail')
        
    elseif(q=='s')
        movefile(dataFilePath, manualPlaceCellStrong)
        movefile(currFilePath, manualPlaceCellStrong)
        disp('strong')
    else
        movefile(dataFilePath, manualPlaceCellMaybe)
        movefile(currFilePath, manualPlaceCellMaybe)
        disp('maybe')
    end
    %data=load(fullfile(dataDir,dataFileName))
    
    %disp('here') 
end

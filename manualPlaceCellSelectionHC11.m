imageDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages';
dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData';

manualPlaceCellFail=fullfile(dataDir,'manuallySortedPlaceCellsAbove5Hz','fail');
manualPlaceCellStrongRightward=fullfile(dataDir,'manuallySortedPlaceCellsAbove5Hz','strongRightward');
manualPlaceCellStrongLeftward=fullfile(dataDir,'manuallySortedPlaceCellsAbove5Hz','strongLeftward');
manualPlaceCellStrongBothWays=fullfile(dataDir,'manuallySortedPlaceCellsAbove5Hz','strongBothWays');

touchDir(manualPlaceCellFail);
touchDir(manualPlaceCellStrongRightward);
touchDir(manualPlaceCellStrongLeftward);
touchDir(manualPlaceCellStrongBothWays)

filePaths=getFilePathsRegex(imageDir,'*tif');

for fi=1:length(filePaths)
    currFilePath=filePaths{fi};
    
    fileName=getFileNameFromPath(currFilePath);
   
    fH = figure(1);
    
    imshow(currFilePath)
    if(fi==1)
     maxFig
    end
    drawnow
    w = waitforbuttonpress;
    q = double(get(gcf,'CurrentCharacter'));
    

    dataFileName=[fileName(1:(end-32)) '.mat'];
    dataFilePath=fullfile(dataDir,dataFileName);
    
    if(q=='d')
        movefile(dataFilePath, manualPlaceCellFail)
        movefile(currFilePath, manualPlaceCellFail)
        disp('fail')
        
    elseif(q=='s')
        movefile(dataFilePath, manualPlaceCellStrongLeftward)
        movefile(currFilePath, manualPlaceCellStrongLeftward)
        disp('strong leftward')
    elseif(q=='a')
        movefile(dataFilePath, manualPlaceCellStrongRightward)
        movefile(currFilePath, manualPlaceCellStrongRightward)
        disp('strong rightward')
    elseif(q=='z')
        movefile(dataFilePath, manualPlaceCellStrongBothWays)
        movefile(currFilePath, manualPlaceCellStrongBothWays)
        disp('strong both ways')
    end
    %data=load(fullfile(dataDir,dataFileName))
    
    %disp('here') 
    
end


dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
imageDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImagesTime';
%imageDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImages/redoCirc';
close all

touchDir(imageDir);
%filePaths=getFilePathsRegex(imageDir,'*tif');
filePaths=getFilePathsRegex(imageDir,'*png');
overwriteFlag=1;

%for fi=1:length(filePaths)
    %for fi=128:length(filePaths)
    for fi=1:1
    currFilePath=filePaths{fi};
    
    fileName=getFileNameFromPath(currFilePath);
    
    currFieldDirection='';
    if(contains(fileName,'leftward'))
        currFieldDirection='leftward';   
    elseif(contains(fileName,'rightward'))
        currFieldDirection='rightward';  
    end

fieldNum=str2num(getSubstrBtwnSubstrs(fileName,'ward_Field','.png'));    
    
    
     firstChar=fileName(1);
    dataFileBaseName=getSubstrBtwnSubstrs(fileName,firstChar,'Info_');
    dataFileName=[firstChar dataFileBaseName 'Info.mat']
    dataFilePath=fullfile(dataDir,dataFileName);
    
    spaceTimeFileName=[firstChar dataFileBaseName sprintf('Info_SpaceTimeResponses_%s.tif',currFieldDirection)];
    
    saveVarName=sprintf('manualTimeField%dEnd%sSec',fieldNum,currFieldDirection);
    data=load(dataFilePath);
    
    if(overwriteFlag==0)
	    if(isfield(data,saveVarName))
		continue
	    end
    end
    
    
    imshow(currFilePath)
    figure
    imshow(spaceTimeFileName)
    showOn2ndMonitor
    %set(gcf,'Position',[-500 500 400 300])
    %if(fi==1)
     %maxFig
    %end
    drawnow
    
    manualTimeFieldEndrightwardSec=[];
    manualTimeFieldEndleftwardSec=[];

    %while(true)
    data
	%data.(saveVarName)
        prompt='Place field/precession end (sec)?';
        manualFieldEndM=input(prompt);
        
	%{
        if(strcmp(currFieldDirection,'rightward'))
             manualTimeFieldEndrightwardSec=[manualTimeFieldEndrightwardSec manualFieldEndM];   
        elseif(strcmp(currFieldDirection,'leftward'))
             manualTimeFieldEndleftwardSec=[manualTimeFieldEndleftwardSec manualFieldEndM];   
        end
	%}
        
    %end
    eval(sprintf('%s = %.2f',saveVarName,manualFieldEndM))
    save(dataFilePath,saveVarName,'-append')
%{
   if(strcmp(currFieldDirection,'rightward'))
     save(dataFilePath,'manualTimeFieldEndrightwardSec','-append')
   elseif(strcmp(currFieldDirection,'leftward'))
     save(dataFilePath,'manualTimeFieldEndleftwardSec','-append') 
   end
 
 %} 
       data=load(dataFilePath)
       fi
end

dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
imageDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImagesTime';
%imageDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImages/redoCirc';
imageDir='/Users/tibinjohn/thetaSeq/code';

close all

touchDir(imageDir);
%filePaths=getFilePathsRegex(imageDir,'*tif');
filePaths=getFilePathsRegex(imageDir,'*_SpaceTimeResponses_*tif');
overwriteFlag=0;

for fi=1:length(filePaths)
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
    
    saveVarName=sprintf('manualTimeField%dEnd%sSec',fieldNum,currFieldDirection);
    data=load(dataFilePath);
    
    if(overwriteFlag==0)
	    if(isfield(data,saveVarName))
		continue
	    end
    end
    
    imshow(currFilePath)
    
    %if(fi==1)
     maxFig
    %end
    drawnow
    
    manualTimeFieldEndrightwardSec=[];
    manualTimeFieldEndleftwardSec=[];

    %while(true)
	data
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

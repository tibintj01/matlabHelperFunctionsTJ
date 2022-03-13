dataDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/manuallySortedPlaceCellsAbove5Hz/strongAllDataFiles';
imageDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImages';
imageDir='/Users/tibinjohn/thetaSeq/code/hc11ThetaCycleData/placeCellQualityImages/singleDirectionGoodPrecessionImages/redoCirc';
close all
%filePaths=getFilePathsRegex(imageDir,'*tif');
filePaths=getFilePathsRegex(imageDir,'*png');
overwriteFlag=1;

for fi=1:length(filePaths)
    currFilePath=filePaths{fi};
    
    fileName=getFileNameFromPath(currFilePath);
    
    currFieldDirection='';
    if(contains(fileName,'leftward'))
        currFieldDirection='leftward';   
    elseif(contains(fileName,'rightward'))
        currFieldDirection='rightward';  
    end
    
    
    
     firstChar=fileName(1);
    dataFileBaseName=getSubstrBtwnSubstrs(fileName,firstChar,'Info_');
    dataFileName=[firstChar dataFileBaseName 'Info.mat']
    dataFilePath=fullfile(dataDir,dataFileName);
    
   
    data=load(dataFilePath);
    
    if(overwriteFlag==0)
    if(isfield(data,sprintf('manualFieldStartsM%s',currFieldDirection)) && isfield(data,sprintf('manualFieldStartsM%s',currFieldDirection)) )
        continue
    end
    end
    
    imshow(currFilePath)
    
    %if(fi==1)
     maxFig
    %end
    drawnow
    
    manualFieldStartsMrightward=[];
    manualFieldEndsMrightward=[];
    
    manualFieldStartsMleftward=[];
    manualFieldEndsMleftward=[];
    

    while(true)
        prompt='Place field/precession start (meters)?';
        manualFieldStartM=input(prompt);

        prompt='Place field/precession end (meters)?';
        manualFieldEndM=input(prompt);
        
        if(strcmp(currFieldDirection,'rightward'))
             manualFieldStartsMrightward=[manualFieldStartsMrightward manualFieldStartM];
             manualFieldEndsMrightward=[manualFieldEndsMrightward manualFieldEndM];   
        elseif(strcmp(currFieldDirection,'leftward'))
             manualFieldStartsMleftward=[manualFieldStartsMleftward manualFieldStartM];
             manualFieldEndsMleftward=[manualFieldEndsMleftward manualFieldEndM];
        end

       
        
        prompt='Enter bounds for another field? (0/1)';
        anotherGoodField=input(prompt);
        
        if(anotherGoodField~=1)
            break
        end
        
    end
     
   if(strcmp(currFieldDirection,'rightward'))
     save(dataFilePath,'manualFieldStartsMrightward','manualFieldEndsMrightward','-append')
   elseif(strcmp(currFieldDirection,'leftward'))
     save(dataFilePath,'manualFieldStartsMleftward','manualFieldEndsMleftward','-append') 
   end
   
       data=load(dataFilePath)
       fi
end

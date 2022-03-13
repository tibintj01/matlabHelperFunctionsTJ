%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot real place cell phase properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
dataDir='/Users/tibinjohn/Downloads/processed_data';
saveDir='exampleHPCCells';
touchDir(saveDir)

filePaths=getFilePathsRegex(dataDir,'*mat');

useExample=0;
%useExample=1;

filterByPhasePrecession=1;
%filterByPhasePrecession=0;

makePlots=1;
makePlots=0;

placeCellPaths=table;

placeCellCount=0;
%for i=1:10
for i=1:length(filePaths)
    tic
    disp(i/length(filePaths))
%for i=281:length(filePaths)
    data=load(filePaths{i});
    
    if(useExample)
        data=load(fullfile(dataDir,'LEM3216_S20190726184722.mat'));
    end
        
    maxDist=max([data.linear_track{1}.lapinfo.laps.pos]);
    posAxis=linspace(0,maxDist,round(maxDist/3));
    meanVelocityPerSession=data.BasicLoco.MeanVelocity;
    numCells=size(data.measures,1);
    numVariables=size(data.measures,2);
    numSessions=size(data.measures,3);
    
    numPanels=3;
    for cellNum=1:numCells
    %for cellNum=3:3
        for sessionID=1:2
            numSessions=1;
          phSlopeTest=abs(max(data.measures(cellNum,14,sessionID)));
         phRSquareTest=abs(max(data.measures(cellNum,15,sessionID)));
         infoContent=max(data.measures(cellNum,1,sessionID));
         
         peakRate=max(data.measures(cellNum,4,sessionID));
         %sparsityTest=max(data.measures(cellNum,3,1:2));
         nFieldsTest=max(data.measures(cellNum,8,sessionID));
         fieldWidth=max(data.measures(cellNum,7,sessionID));
         lapPermStability=data.measures(cellNum,55,sessionID);
         nSpikes=data.measures(cellNum,9,sessionID);
         
         if(nFieldsTest==0)
             disp('ZEROO')
             %break
         end
         %disp(infoContent)
         if(filterByPhasePrecession)
             %if((phSlopeTest)<4 || (phRSquareTest)<0.3 || isnan(phSlopeTest) ||...
             %        infoContent<0.2 || peakRate<5 || nFieldsTest==0)
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %criteria for place cells from Harvey et al., 2020
             %except criterion for "15 laps with consistent behavior"
             %(lapPermStability>0.3 also included) make 0 to exclude this
             %criterion
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             min_lapPermStability=0.3;
             %min_lapPermStability=0;
             %if(peakRate<1 || fieldWidth <9 || fieldWidth>78 || lapPermStability<min_lapPermStability ...
             if(peakRate<10 || fieldWidth <9 || fieldWidth>78 || lapPermStability<min_lapPermStability ...
                     || nSpikes<100 || infoContent<0.15)
                    %%|| nSpikes<100 || infoContent<0.15)
                 %close all
                 continue
             end

	     if(isnan(peakRate) || isnan(fieldWidth) || isnan(lapPermStability) || isnan(nSpikes) || isnan(infoContent))
		continue
	     end 
         end
         %record the path leading to place cell criterion being passed
         fileNameRoot=getFileNameFromPath(filePaths{i});
        fileNameRoot=fileNameRoot(1:end-4);
        placeCellCount=placeCellCount+1;
        
        %data.measures(cellNum,50,sessionID); %spike amplitude
        
        tableRow=table({fileNameRoot},cellNum,sessionID);
        placeCellPaths=[placeCellPaths; tableRow];
        
 
        if(makePlots)
            plotPlaceCellPropertiesClark
        end
        
        end
    end
    if(useExample)
        break
    end
    toc
end
placeCellPaths.Properties.VariableNames={'FileName','cellNum','sessionNum'};
%save(fullfile(saveDir,'placeCellFilePaths.mat'),'placeCellPaths','placeCellCount','min_lapPermStability')
save(fullfile(saveDir,'placeCellFilePathsMinPeakRate10.mat'),'placeCellPaths','placeCellCount','min_lapPermStability')
findNumPlaceCells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data from simulation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010/fluxRawData/stabilityData';
inputDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010/';

close all

minPerturbTime=1;
maxPerturbTime=150;
nsyn=7
directionStr='distToProx';
%directionStr='proxToDist';
%precessionType='log';
precessionType='linear';

numStatesPerDim=500
minV=-100;
maxV=-40;

numBatches=20;
runsPerBatch=numStatesPerDim*numStatesPerDim/numBatches;

totalCount=numStatesPerDim*numStatesPerDim;

saveDt=0.5;
%slopeEstDeltaT=10; %msec
slopeEstDeltaT=5; %msec
slopeEstDeltaT=3; %msec
slopeEstDeltaIdx=round(slopeEstDeltaT/saveDt);
count=0;

sampleStep=10;
sampleStep=15;
x=NaN(round(totalCount/sampleStep),1);
y=NaN(round(totalCount/sampleStep),1);
u=NaN(round(totalCount/sampleStep),1);
v=NaN(round(totalCount/sampleStep),1);

for batchNum=1:(numBatches-1)
    currInputFileName=sprintf('perturbationConfigsSection%d.dat',batchNum)
    currFileName=sprintf('tempV_CenterSyn0_section%d_%sPrecession.dat',batchNum,precessionType);
   
    currInput=readtable(fullfile(inputDir,currInputFileName));
    data=readtable(fullfile(dataDir,currFileName));
     if(isempty(currInput))
         continue
     end
    
    for currRun=1:sampleStep:runsPerBatch
            count=count+1;
                
        currV=data{1:(end-1),currRun};
       
        
        currInitialPtT=currInput{currRun,1};
        currInitialPtV=currInput{currRun,2};
        
        
        
            
        currPerturbTime=minPerturbTime+(maxPerturbTime-minPerturbTime)*currInput{currRun,1};
        
       
        currPerturbIdx=round(currPerturbTime/saveDt)+3;
        %{
        currPerturbIdx=round(currPerturbTime/saveDt);
        baselineVal=currV(currPerturbIdx-1)
        
        searchStartIdx=max(1,currPerturbIdx-slopeEstDeltaIdx/5);
        searchEndIdx=min(length(currV),currPerturbIdx+slopeEstDeltaIdx/5);
        [~,maxIdxInc]=max(abs(currV(searchStartIdx:searchEndIdx)-baselineVal));
        searchWindowLength=length(searchStartIdx:searchEndIdx);
       currPerturbIdx=currPerturbIdx+round(maxIdxInc-searchWindowLength/2)+1;
       %}
        %currPerturbIdx=min(length(currV),currPerturbIdx);
        
        
        %fds
          
%figure
%plot(currV)
        currPerturbFinalIdx=currPerturbIdx+slopeEstDeltaIdx;
        currPerturbFinalIdx=min(length(currV),currPerturbFinalIdx);
       
        initialV=currV(currPerturbIdx);
        finalV=currV(currPerturbFinalIdx);
        
        currDeltaV=finalV-initialV;
        currDeltaT=slopeEstDeltaT;
        
        x(count)=currPerturbTime;
        y(count)=initialV;
        u(count)=currDeltaT;
        v(count)=currDeltaV;

        %figure(1)
        %plot(currV)
        %hold on
    end
    count
    figure
    %quiver(x,y,u,v,3)
     quiver(x,y,u,v,2.5)
    xlabel('Time in theta cycle (msec)')
    ylabel('V (mV)')
    maxFig
    setFigFontTo(18)
    drawnow
    
    saveas(gcf,'sequenceStabilityVectorField.tif')
    close all
end

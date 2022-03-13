%numStatesPerDim=10;
numStatesTimeDim=500;
numStatesVDim=500;
%numStatesTimeDim=50;
%numStatesVDim=50;
numPerturbations=numStatesTimeDim*numStatesVDim;

numBatches=20;

countPerBatch=floor(numPerturbations/numBatches);

perturbationIDs=1:numPerturbations;

saveDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010';

perturbationPts=zeros(numPerturbations,2);

count=1;
for x1=1:numStatesTimeDim
    for x2=1:numStatesVDim
        currStateVec=[x1 x2];
        
        perturbationPts(count,:)=(currStateVec-0.5)/numStatesTimeDim;
        
        if(mod(count,countPerBatch)==0 || count==numPerturbations)
             if(count==numPerturbations)
                  startIdx=countPerBatch*numBatches+1;
             else
                startIdx=count-countPerBatch+1;
             end
            dlmwrite(fullfile(saveDir,sprintf('perturbationConfigsSection%d.dat',ceil(count/countPerBatch))),perturbationPts(startIdx:count,:),'delimiter',' ')
        end
        count=count+1;    
    end
end
perturbationPts=perturbationPts(1:(count-1),:);

%dlmwrite(fullfile(saveDir,'statesInSphere.dat'),spherePts,'delimiter',' ')

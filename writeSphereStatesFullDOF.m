%numStatesPerDim=10;
numStatesPerDim=8;
dof=7;
numSphereStates=numStatesPerDim^dof;

numBatches=20;

countPerBatch=floor(numSphereStates/numBatches);

sphereStateIDs=1:numSphereStates;

saveDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010';

spherePts=zeros(numSphereStates,dof);

count=1;
for x1=1:numStatesPerDim
%x1=1
x1
    for x2=1:numStatesPerDim
        for x3=1:numStatesPerDim
        	for x4=1:numStatesPerDim
        		for x5=1:numStatesPerDim
        			for x6=1:numStatesPerDim
        				for x7=1:numStatesPerDim
                            currStateVec=[x1 x2 x3 x4 x5 x6 x7];
                            currStateVec(currStateVec==1)=-0.5; %first state is -1 == not active
		     				spherePts(count,:)=(currStateVec-0.5)/numStatesPerDim ;
                            if(mod(count,countPerBatch)==0 || count==numSphereStates)
                                 if(count==numSphereStates)
                                      startIdx=countPerBatch*numBatches+1;
                                 else
                                    startIdx=count-countPerBatch+1;
                                 end
                                dlmwrite(fullfile(saveDir,sprintf('statesInSphereSection%d.dat',ceil(count/countPerBatch))),spherePts(startIdx:count,:),'delimiter',' ')
                            end
		     				count=count+1;
                        end
                    end
                end
            end
        end
    end
end
spherePts=spherePts(1:(count-1),:);

%dlmwrite(fullfile(saveDir,'statesInSphere.dat'),spherePts,'delimiter',' ')

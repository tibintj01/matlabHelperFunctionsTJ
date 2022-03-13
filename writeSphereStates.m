%numStatesPerDim=10;
numStatesPerDim=10;
numSphereStates=numStatesPerDim^3;


sphereStateIDs=1:numSphereStates;

saveDir='/Users/tibinjohn/thetaSeq/code/hocFiles/Branco_2010';

spherePts=zeros(numSphereStates,3);

count=1;
for xi=1:numStatesPerDim
    for yi=1:numStatesPerDim
        for zi=1:numStatesPerDim
            %if(abs(xi-yi) <= 5 || abs(xi-zi) <= 5 ||abs(yi-zi) <= 5)
            %    continue
            %end
            %{
            if(xi>numStatesPerDim/3)
                continue
            end
            
             if(yi>2*numStatesPerDim/3 || yi<numStatesPerDim/3)
                continue
             end
              if(zi<2*numStatesPerDim/3)
                continue
            end
            %}
            
             spherePts(count,:)=([xi yi zi]-0.5)/numStatesPerDim ;
             count=count+1;
        end
    end
end
spherePts=spherePts(1:(count-1),:);

dlmwrite(fullfile(saveDir,'statesInSphere.dat'),spherePts,'delimiter',' ')
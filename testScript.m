%[W,H]=seqNMFaccessRoutine(xValues,A,columnDescrip,xUnitStr,iterationUnitStr,yUnitStr)

%trajectory is distance from each component at each time
numFactorsUsed=17;

assemblySimilarityOverTime=NaN(numFactorsUsed,cycleID);
for i=1:cycleID
        thisCycleActivity=A(:,i)'; %1xN
        assemblyPatternMatrix=W(:,1:numFactorsUsed);%Nx18

        thisCyclePosition=thisCycleActivity*assemblyPatternMatrix;
        assemblySimilarityOverTime(:,i)=thisCyclePosition(:);
end
figure
imagesc(assemblySimilarityOverTime)
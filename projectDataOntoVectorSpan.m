function [projectedData] = projectDataOntoVectorSpan(dataMatrix,desiredVector)
%projects data matrix (each datapt is a col) onto span of given vector

    %1xN * NxM= 1xM (# of datapts)
    projectedData=desiredVector(:)'*dataMatrix; 
end


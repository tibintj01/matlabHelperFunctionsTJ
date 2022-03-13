function [activationVec] = listToActivationVec(listIn,maxVal)
    activationVec=zeros(maxVal,1);
    
    for i=1:length(listIn)
        activationVec(listIn(i))=i;
    end
    
end


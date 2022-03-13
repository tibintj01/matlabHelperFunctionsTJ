function [H] = getEntropy(p)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    H=0;
    for i=1:length(p)
        if(p(i)>0)
            H=H-p(i)*log(p(i))/log(2);
        end
    end
end


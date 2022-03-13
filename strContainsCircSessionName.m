function [isCircMaze] = strContainsCircSessionName(currFileName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
circularSessionNames={'Achilles_11012013','Gatsby_08282013','Cicero_09102014'};
    isCircMaze=0;
    if((contains(currFileName,circularSessionNames{1}) || contains(currFileName,circularSessionNames{2}) || contains(currFileName,circularSessionNames{3})))
        isCircMaze=1;
    end
end


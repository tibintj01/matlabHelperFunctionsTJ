function [currSessionNum,currSessName]=getSessionNum(fileBaseName)
%UNTITLED Summary of this function goes here
%sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09102014','Cicero_09172014','Gatsby_08282013'};
sessionNames={'Achilles_10252013','Achilles_11012013','Buddy_06272013','Cicero_09012014','Cicero_09102014','Cicero_09172014','Gatsby_08282013'};


currSessName=NaN;
currSessionNum=0;
for i=1:length(sessionNames)
    if(contains(fileBaseName,sessionNames{i}))
        currSessionNum=i;
        currSessName=sessionNames{i};
        break
    end
end



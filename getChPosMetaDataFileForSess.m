function [chPosMetaDataThisSessPath]=getChPosMetaDataFileForSess(channelPosMetaDataDir,currSessionName)
%find ch pos file corresponding to this session name

chPosMetaDataThisSessPath=NaN;

if(contains(currSessionName,'.'))
    animalName=getStrUpTo(currSessionName,'.');
else
    animalName=currSessionName;
end



candidateFileNames=getFilePathsRegex(channelPosMetaDataDir,sprintf('%s*MaxPosi.txt',animalName));

for fi=1:length(candidateFileNames)
    currCandidateFilePath=candidateFileNames{fi};
    currCandidateFile=getFileNameFromPath(currCandidateFilePath);
    
    %if(contains(currSessionName,'vvp') || contains(currSessionName,'gor'))
    if(~contains(currSessionName,'ec'))
        fileNameRoot=getStrUpTo(currCandidateFile,'.')
        
        if(contains(currSessionName,fileNameRoot))
             chPosMetaDataThisSessPath=currCandidateFilePath;
            break
        end
        
    else
        currSessionNum=str2num(getSubstrBtwnSubstrs(currSessionName,'.',''));
        currCandidateFileBaseName=getStrUpTo(currCandidateFile,'M');
        fileNameRoot=getSubstrBtwnSubstrs(currCandidateFileBaseName,'.','.');
        
        startSesNum=str2num(fileNameRoot(1:3));
        endSesNum=str2num(fileNameRoot(5:7));
        
        if(currSessionNum>=startSesNum && currSessionNum<=endSesNum)
            
            chPosMetaDataThisSessPath=currCandidateFilePath;
            break
        end
        
    end
    
    
    
end



function [fileName] = getFileName(fileNamePath)
        pathSplit=strsplit(fileNamePath,'/');
        if(length(pathSplit)==1)
                pathSplit=strsplit(fileNamePath,'\');
        end
        fileName=pathSplit{end};

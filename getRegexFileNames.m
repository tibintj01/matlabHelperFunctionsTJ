function [fileNames]=getRegexFilePaths(dirName,regexStr)
                dirInfo=dir(fullfile(dirName,regexStr));
                fileNames={dirInfo(:).name};

function [fileName]=getRegexFilePaths(dirName,regexStr)
                dirInfo=dir(fullfile(dirName,regexStr));
                fileNames={dirInfo(:).name};

		if(~isempty(fileNames))
			fileName=fileNames{1};
		else
			fileName='RatingTooLow';
		end

function [fileName] = getDirName(fileNamePath)
	pathSplit=strsplit(fileNamePath,'/');
	if(length(pathSplit)==1)
		pathSplit=strsplit(fileNamePath,'\');	
	end
	fileName=fullfile(pathSplit{1:(end-1)});
	if(fileNamePath(1)=='/')
		fileName=['/' fileName];
	end

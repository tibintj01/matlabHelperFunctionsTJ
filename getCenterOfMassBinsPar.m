function [] = getCenterOfMassBinsPar(startIdx,endIdx) 
tic
if(isstr(startIdx))
	startIdx=str2num(startIdx);
end
if(isstr(endIdx))
	endIdx=str2num(endIdx);
end
filePaths=getRegexFilePaths('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties','*.mat');

for i=startIdx:endIdx
	i
	[troughCMBin,peakCMBin]=getCentersOfMass(filePaths{i},0.001,0.001,pwd);
	fileName=getFileNameFromPath(filePaths{i});
	save(sprintf('centerOfMassBins%s',fileName),'troughCMBin','peakCMBin')

end
disp('Single run time:')
toc



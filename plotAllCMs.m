filePaths=getRegexFilePaths('/nfs/turbo/lsa-ojahmed/compiledDir_01-Nov-2017_23-47-00_getCellWaveformPropsOmar_Tibin/seizureCellWaveformProperties','*.mat')

for i=1:length(filePaths)
	getCentersOfMass(filePaths{i},0.001,0.001,pwd)
end



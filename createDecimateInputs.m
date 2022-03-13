decimateDir='/nfs/turbo/lsa-ojahmed/FOR_TIBIN_otherPts/MG29/NS5/';

decimateFileNames=getRegexFileNames(decimateDir,'*ns5');

for i=1:length(decimateFileNames)
	inputCellArray{i,1}=decimateDir;
	inputCellArray{i,2}=decimateFileNames{i};
end

save('mg29decimateInputs.mat','inputCellArray')

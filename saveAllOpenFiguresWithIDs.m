function [] = saveAllOpenFiguresWithIDs(saveDir,figNameList)
%FolderName = pwd;   % Your destination folder
FolderName = saveDir;   % Your destination folder
touchDir(FolderName)

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  %FigName   = get(FigHandle, 'Name');
   %FigName=sprintf('Fig%d.tif',iFig);
   FigName=sprintf(sprintf('%s.tif',figNameList{iFig}));
  saveas(FigHandle, fullfile(FolderName, FigName));
end

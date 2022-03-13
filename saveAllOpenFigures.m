function [] = saveAllOpenFigures(filePath)

 figHandles = findall(0,'Type','figure'); 
 
 if(isempty(figHandles))
     disp('no open figures')
 else
     % Save first figure
 export_fig(filePath, '-pdf', figHandles(1))
 
 % Loop through figures 2:end
 for i = 2:numel(figHandles)
     export_fig(filePath, '-pdf', figHandles(i), '-append')
 end
 
 end

function [row,col]= getSpatialRowColFromFileName(fileName,ch)

	fileName=getFileNameFromPath(fileName);
        subject=getSubject(fileName);
	
	[arrayMap,electrodeXY,electrodeImp]=neuroportArrayData(subject);

	%electrodeXY
	%row=electrodeXY(ch,1);
	%col=electrodeXY(ch,2);
	
	col=electrodeXY(ch,1);
	row=(10-electrodeXY(ch,2))+1;
end

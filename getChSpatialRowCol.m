function [row,col]= getChSpatialRowCol(ptDir,ch)

	[arrayMap,electrodeXY,electrodeImp]=neuroportArrayData(ptDir);

	
	%electrodeXY
	%row=electrodeXY(ch,1);
	%col=electrodeXY(ch,2);
	
	col=electrodeXY(ch,1);
	row=(10-electrodeXY(ch,2))+1;
end

function [row,col,cellIdx]= getSpatialRowColFromFileName(fileName)
	fileName=getFileNameFromPath(fileName);
	subject=getSubject(fileName);
	cellIdx='x';
	%chCell=getChCell(fileName);
	%ch=str2num(getStrUpTo(chCell,'_'));
	chCell=(getStrUpTo(fileName,'_'));
	if(length(chCell)==2)
		ch=str2num(chCell(1));
		%cellIdx=letter2num(chCell(2));
	elseif(length(chCell)==3)
		ch=str2num(chCell(1:2));
		%cellIdx=letter2num(chCell(3));
	end	

	ch
	
	[arrayMap,electrodeXY,electrodeImp]=neuroportArrayData(subject);

	electrodeXY
	%row=electrodeXY(ch,1);
	%col=electrodeXY(ch,2);
	
	col=electrodeXY(ch,1);
	row=(10-electrodeXY(ch,2))+1;
end

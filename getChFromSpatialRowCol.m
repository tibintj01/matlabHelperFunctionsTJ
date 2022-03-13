function [chNum]=getChFromSpatialRowCol(cellPropFile,row,col)
                fileName=getFileNameFromPath(cellPropFile);
                subject=getSubject(fileName);

                [arrayMap,electrodeXY,electrodeImp]=neuroportArrayData(subject);

		chNum=-1;
		for ch=1:96
			if(electrodeXY(ch,1)==col && ((10-electrodeXY(ch,2))+1 ==row))
				chNum=ch;
			end
		end

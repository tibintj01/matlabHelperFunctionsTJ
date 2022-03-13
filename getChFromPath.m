function [ch] = getChFromPath(filePath)

	fileName=getFileNameFromPath(filePath);

	chCell=(getStrUpTo(fileName,'_'));
        if(length(chCell)==2)
                ch=str2num(chCell(1));
                %cellIdx=letter2num(chCell(2));
        elseif(length(chCell)==3)
                ch=str2num(chCell(1:2));
                %cellIdx=letter2num(chCell(3));
        end	

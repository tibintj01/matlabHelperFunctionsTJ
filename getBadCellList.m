function [badCells]=getBadCellList(cellPropFileName)
	%bad = ranked less than 2
	badCells={''};
	if(length(findstr(cellPropFileName,'MG49-seizure45'))>0)
                %badCells={'3a','8a','12c','16a','22a','22b','27b','31c','33a','37a','47c','48a','49a'0
                badCells={'12c','22a','37a','47c','62b','66a','74a','76a','78a'};
        end
        if(length(findstr(cellPropFileName,'MG49-seizure43'))>0)
        end
        if(length(findstr(cellPropFileName,'MG49_seizure36'))>0)
        end
        if(length(findstr(cellPropFileName,'MG63_seizure1-4'))>0)
        end
        if(length(findstr(cellPropFileName,'BW9_seizure1-3'))>0)
        end
        if(length(findstr(cellPropFileName,'RIHE1'))>0)
        end	

function [xgram]=tryAll2DShifts(img1,img2,maxShiftCol,maxShiftRow)
%gets limited range crosscorrelogram

 
    colShifts=(-maxShiftCol):maxShiftCol;
    rowShifts=(-maxShiftRow):maxShiftRow;
       colShifts=(0):maxShiftCol;
    rowShifts=(0):maxShiftRow;
    
       xgram=NaN(length(rowShifts),length(colShifts));
       
    for ci=1:length(colShifts)
        ci
        for ri=1:length(rowShifts)
            currColShift=colShifts(ci);
            currRowShift=rowShifts(ri);
            
            shiftedImg2=imtranslate(img2,[currRowShift currColShift]);
            
            currCorrVal=corrcoef([img1(:) shiftedImg2(:)]);
            currCorrVal=currCorrVal(1,2);
            
            xgram(ri,ci)=currCorrVal;
        end
    end




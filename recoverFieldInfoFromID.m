function [ui,di,dfi]=recoverFieldInfoFromID(currFieldID)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UNIQUE FIELD ID CODE: unitIDinSession*1000+di*100+directionalfieldNum
%this code recovers the original variables that went into this code^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fieldIDStr=num2str(currFieldID);

%places=[1 10 100 1000 10000];
places=[10000 1000 100 10 1];
digits=NaN(size(places));

for digitI = 1:length(places)
    d=places(digitI);
    digits(digitI)=round( (mod(currFieldID,10*d)-mod(currFieldID,d))/d );
end


dfi=digits(end);
di=digits(3);

if(digits(1)>0)
    ui=str2num(sprintf('%d%d',digits(1),digits(2)));
else
    ui=digits(2);
end










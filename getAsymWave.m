function [f] = getAsymWave(t,troughFrac,amp)
%UNTITLED2 Summary of this function goes here
%maxT=5;
%troughFrac=0.7;
if(~exist('amp','var'))
amp=2;
end

%stepSize=0.01;
%t=0:stepSize:maxT;

f=NaN(size(t));
expDischargeOfInh=0;

for i=1:length(t)
    currT=mod(t(i),1); 
    if currT <= troughFrac
        if(expDischargeOfInh)
            f(i) =  getLogTime(1-currT,1.4);
        else
            f(i) = -1 / troughFrac * currT + 1; 
        end
    else  
        f(i) = 1 / (1-troughFrac) * (currT - troughFrac);
    end
end

%f=f-median(f);



%upSampStep=stepSize/10;

%f=interp1(t,f,0:upSampStep:maxT);

f=f.*amp;

disp('')


%figure; plot(t,f,'k-','LineWidth',4)

%hold on



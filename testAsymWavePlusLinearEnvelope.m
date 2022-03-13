close all; clear all; clc

maxT=10;
maxT=20;
maxT=50;
maxT=20;
maxT=15;
maxT=10;

troughFrac=0.7;
%troughFrac=0.9;

useLinExc=0;
%useLinExc=1;

stepSize=0.01;
t=0:stepSize:maxT;

f=NaN(size(t));

for i=1:length(t)
    currT=mod(t(i),1); 
    if currT <= troughFrac
        f(i) = -1 / troughFrac * currT + 1;  
    else  
        f(i) = 1 / (1-troughFrac) * (currT - troughFrac);
    end
end

f=f-median(f);
f=f*2;


upSampStep=stepSize/10;
f=interp1(t,f,0:upSampStep:maxT);
f=scaledata(f,0.3,1);
%f=scaledata(f,0,0.5);

%excIn=scaledata(excIn,0,1);
%excIn=excIn-0.5;
t=0:upSampStep:maxT;
linEnv=t/maxT;
offset=pi/2;
offset=3*pi/4;
offset=pi/2;
offset=pi;
%offset=0;
offset=pi/4;
offset=pi/2;
%offset=-pi;
%offset=-pi/1.5;
if(useLinExc)
    excIn=linEnv;
else
    %excIn=2*linEnv.*(sin(2*pi*t+offset)+1)/2;
    excIn=2*getAsymWave(t+0.25,troughFrac,linEnv);
     %excIn=2*(linEnv.^2).*(sin(2*pi*t+offset)+1)/2;

end


excInSlope=[diff(excIn) NaN];
inhSlope=[diff(f) NaN];
%excIn=1*linEnv;



%ylim([0 1])



tolLevel=upSampStep*3;
crossIdxes=find(abs(f-excIn)<tolLevel & (excInSlope>inhSlope)>0);
figure; 
subplot(2,1,1)
plot(t,f,'k-','LineWidth',1)
hold on

plot(t,excIn,'b-','LineWidth',1)
hold on
plot(t(crossIdxes),f(crossIdxes),'r.','MarkerSize',40)

legend('asym inh','rhythmic exc','exc surpasses inh','Location','northwest')
legend boxoff
box off
xlabel('Time')
ylabel('Amplitude')

ylim([0 2])

subplot(2,1,2)
plot(t(crossIdxes),mod(t(crossIdxes),1),'ro')
ylim([0 1])

setFigFontTo(18)
maxFig

if(useLinExc)

saveas(gcf,'asymThetaPhasePrecessionModelLinExc.png')
else
    saveas(gcf,'asymThetaPhasePrecessionModelRhythmExc.png')
end
%xlim([25 Inf])

close all; clear all; clc

maxT=10;
maxT=20;
maxT=50;
maxT=20;
maxT=15;
troughFrac=0.9;
troughFrac=1;



stepSize=0.01;
t=0:stepSize:maxT;

f=NaN(size(t));

useLinExc=1;
%useLinExc=0;

expDischargeOfInh=1;

expDischargeOfInh=1;

for i=1:length(t)
    currT=mod(t(i),1); 
    if currT <= troughFrac
        if(expDischargeOfInh)
            %f(i) =  getLogTime(1-currT,1.4);
             %f(i) =  getLogTime(1-currT,1.2);
            f(i) =  getLogTime(1-currT,1.7);
        else
            f(i) = -1 / troughFrac * currT + 1; 
        end
    else  
        f(i) = 1 / (1-troughFrac) * (currT - troughFrac);
    end
end

f=f-median(f);
f=f*2;


upSampStep=stepSize/10;

f=interp1(t,f,0:upSampStep:maxT);


f=scaledata(f,0,1);
%f=scaledata(f,0,0.5);

%excIn=scaledata(excIn,0,1);
%excIn=excIn-0.5;
t=0:upSampStep:maxT;
linEnv=t/maxT;
offset=pi/2;
offset=3*pi/4;

offset=pi;
%offset=0;
offset=pi/4;
offset=0;
offset=-pi/2;
if(useLinExc)
    amp=100;
    amp=20;
    growthRate=5;
    amp=20;
    growthRate=5*.8;
    
    %excIn=(exp(5*linEnv)-1)/100;
    %excIn=(exp(growthRate*linEnv)-1)/amp;
    excIn=linEnv;
else
    excIn=2*linEnv.*(sin(2*pi*t+offset)+1)/2;
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




if(useLinExc)
    legend('concave inh','linear exc','exc overtakes inh','Location','northwest')
else
legend('concave inh','rhythmic exc','exc overtakes inh','Location','northwest')
end
legend boxoff
box off
xlabel('Time (number of cycles elapsed)')
ylabel('Amplitude')

ylim([0 2])

subplot(2,1,2)
plot(t(crossIdxes),mod(t(crossIdxes),1),'r.','MarkerSize',40)
ylim([0 1])

xlabel('Time (number of cycles elapsed)')
ylabel('Theta phase of crossing')

setFigFontTo(18)


maxFigHalfWidth

if(useLinExc)

saveas(gcf,'asymThetaPhasePrecessionModelLinExc.png')
else
    saveas(gcf,'asymThetaPhasePrecessionModelRhythmExc.png')
end
%xlim([25 Inf])

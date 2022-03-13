
close all
x=-0.5:0.001:0.5;

xLog=0:0.001:1;
%xLog=0:0.001:2;

%xExp=0:0.001:10;
xExp=0:0.001:1;



negIdxes=x<0;

%y=(log(x)/log(2)+7)*360/7;

%y=(log(x)/log(2)+7)/7;

%y=logTime(x-0.025);
%y=logTime(x);
%y=logTime(x-.1,);
%y=getLogTime(x+0.05,sqrt(2));
%yLog=getLogTime(1-xLog*.925+0.01,sqrt(2));

base=sqrt(2);
yLog=getLogTime(1-xLog,base);
%yLog=getLogTime(xLog,sqrt(2));

yLog=scaledata(yLog,0,1);

yExp=1-(base.^(xExp*7));
yExp=scaledata(yExp,0,1);

%y=getLogTime(x,sqrt(2));
%y=log(abs(x/10));
%y(negIdxes)=-1*y(negIdxes)-18;

figure
subplot(2,2,1)
plot(xLog,yLog,'k-','LineWidth',3)
box off
%ylim([0 1])
%xlim([0 1])

xlabel('x (a.u.)')
ylabel('y=log(1-x) (a.u.)')

title('log(1-x)')

subplot(2,2,2)
%plot(yExp,xExp,'k-','LineWidth',3)
plot(xExp,yExp,'k-','LineWidth',3)
box off
%ylim([0 360])
%xlim([0 1])

title('1-exp(x)')
xlabel('x (a.u.)')
ylabel('y=1-exp(x) (a.u.)')

ylim([0 1])
xlim([0 1])

subplot(2,2,3)
plot(yLog,([diff(yLog) 0]),'k-','LineWidth',3)
%plot(xLog,1./([diff(yLog) 0]),'k-','LineWidth',3)

box off
hold on

%yyaxis right
%plot(yLog,1./(1-yLog),'k-','LineWidth',3)
%plot(xLog,([diff(yLog) 0]),'k-','LineWidth',3)


xlim([0 1])
xlabel('y=log(1-x) (a.u.)')
ylabel('diff(y) (a.u.)')
ylim([-Inf Inf])

subplot(2,2,4)
plot(yExp,[diff(yExp) 0],'k-','LineWidth',3)

axis tight
box off
xlabel('y=1-exp(x) (a.u.)')
ylabel('diff(y) (a.u.)')

setFigFontTo(18)



%xlim([0 1])


%negsemilogx(x,y)

%symlog(gca,'x')
%set(gca,'xscale','log')
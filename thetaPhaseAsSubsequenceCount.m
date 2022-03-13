close all
%x=0:50;
x=0:100;
x=0:200;
x=0:128;
maxPow=11.5;
maxPow=10;
maxPow=8;
maxPow=7;
maxPow=6.5;
maxPow=6.75;
maxPow=6.85;
maxPow=7;
x=0:0.1:(2^maxPow);

figure
subplot(2,2,1)
logH=plot(x,log2(x)*(7/maxPow),'k','LineWidth',4)

hold on
linH=plot(x,x/(x(end)/7),'b','LineWidth',4)

for i=1:7
    hold on
   plot([x(1)+0.0001 x(end)], [i i],'k--')
end


xlim([1 x(end)])
ylim([0 7])
pbaspect([1.29 1 1]) %same as Mehta et al., 2002
legend([logH linH],{'log_2(x)','linear(x)'},'Location','Best')
title({'Distance (~amount of information contributed to next sequence), linear scale'})

xlabel('Distance to cue (a.u.)')
ylabel('Theta Phase Bin')

subplot(2,2,3)

logH=plot(x,log2(x)*(7/maxPow),'k','LineWidth',4)
title({'Distance (~amount of information contributed to next sequence), log scale','Info-content equalizing representation'})
hold on
linH=plot(x,x/(x(end)/7),'b','LineWidth',4)

for i=1:7
    hold on
   plot([x(1)+0.0001 x(end)], [i i],'k--')
end


xlim([1 x(end)])
ylim([0 7])
pbaspect([1.29 1 1]) %same as Mehta et al., 2002
set(gca,'xscale','log')

legend([logH linH],{'log_2(x)','linear(x)'},'Location','Best')
%legend([logH linH],{'log(x)','linear(x)'})

xlabel('Distance to cue (a.u.)')
ylabel('Theta Phase Bin')

subplot(2,2,[2 4])

binCollectionBias=NaN(7,1);

for i=1:7
   currBinCount=length(find( log2(x)*(7/maxPow) > (i-1) & log2(x)*(7/maxPow) < (i)));
   binCollectionBias(i)=currBinCount/length(x);
end


bar(1:7,binCollectionBias)
xlabel('Theta phase bin')
ylabel('Fraction of curve in bin')

title('Exponential decrease in spike-grouping across theta cycle')

uberTitle({'Logarithmic precession as correcting for info content disparity of place cell spikes', ' across theta cycle about upcoming theta sequence message'},20)



sigma=0.05;
sigma=0.1;
%sigma=0.025;
sigma=0.2;

x0=sigma:0.001:(1-sigma);

mappedIntervalLength=NaN(size(x0));
for i=1:length(x0)
    startInterval=x0(i)-sigma;
    endInterval=x0(i)+sigma;
    
    mappedIntervalLength(i)=logTime(endInterval)-logTime(startInterval);
end

figure; plot(x0,sqrt(mappedIntervalLength),'LineWidth',5)

%set(gca,'yscale','log')
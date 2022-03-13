clear all; close all
x=0:0.01:10; 
p=getGaussCurve(x,5,1);

%p=linspace(1,1,length(x));

p=p/sum(p); %convert to probability measure by sum normalization

figure; plot(x,p,'LineWidth',3)
yyaxis right; plot(x,-log2(p),'LineWidth',3)

infoQuant=0;
for i=1:length(x)
   infoQuant=infoQuant-(p(i)*log2(p(i)));
end

hold on

plot([x(1) x(end)],[infoQuant infoQuant],'k','LineWidth',3)
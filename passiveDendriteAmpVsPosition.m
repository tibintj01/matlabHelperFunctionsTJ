figure

%x=0:0.0001:4;
maxDistFromSynapse=6;
x=-maxDistFromSynapse:0.0001:maxDistFromSynapse;

s1=exp(-(abs(x)));
s2=exp(-(abs(x-0.5)));
s3=exp(-(abs(x-1)));
s4=exp(-(abs(x-1.5)));
s5=exp(-(abs(x-2)));
s6=exp(-(abs(x-2.5)));
s7=exp(-(abs(x-3)));

plot(x,s1,'k')
hold on
plot(x,s2,'b')
plot(x,s3,'r')
plot(x,s4)
plot(x,s5)
plot(x,s6)
plot(x,s7)

plot([4 4],[0 1],'k--')


xlabel('distance along dendrite from synapse (space constants \lambda)')
ylabel('PSP magnitude')
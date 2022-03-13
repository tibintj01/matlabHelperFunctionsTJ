x=linspace(0,1,100);
t=linspace(0,1,100);

[X,T]=meshgrid(x,t);
modelP=(2-(X+T))/2;

figure;

cutOff=0.2;
cutOff=0.3;
modelP(T-X>cutOff)=NaN;
modelP(X-T>cutOff)=NaN;
%modelP(T-X>1-cutOff)=NaN;
%modelP(X>1-cutOff & T<cutOff)=NaN;

fH=figure
subplot(2,2,2)
omarPcolor(x,t,modelP',fH)
cb=colorbar
ylabel(cb,'P(x,t)')
colormap(gca,jet)
xlabel('x')
ylabel('t')
title('diagonal isophase line model')
daspect([1 1 1])
subplot(2,2,1)
plot(x,nanmean(modelP,1))
xlabel('x')
ylabel('P_x')
title('diagonal model x projection')
daspect([1 1 1])

subplot(2,2,4)
plot(t,nanmean(modelP,2))
xlabel('x')
ylabel('P_t')
title('diagonal model t projection')
setFigFontTo(18)
daspect([1 1 1])
numPts=500;
x=linspace(0,1,numPts);
t=linspace(0,1,numPts);

[X,T]=meshgrid(x,t);

P=max(x)^2-X.*T;

figure
subplot(1,2,1)
s=surf(X,T,P)
daspect([1 1 1])
colormap(jet)
colorbar
set(s, 'EdgeColor', 'none')

subplot(1,2,2)
contour(X,T,P,30)
daspect([1 1 1])
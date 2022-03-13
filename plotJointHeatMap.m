function []=plotJointHeatMap(x,y,binWidth,fH)

%assumes y is circ variable

[ rho,p,s,b ] = kempter_lincirc( x,y*2*pi);

nBinsX=round(range(x/binWidth));
nBinsY=round(range(y/binWidth));

    [N,c]=hist3([x(:) y(:)],[nBinsX nBinsY]);
     omarPcolor(c{1},c{2},N',fH)
          cb=colorbar
         ylabel(cb,'Joint probability')
         colormap(gca,jet)
         title(sprintf('circ R= %.2f, p=%.3f',rho,p))
         %{
         xlabel('Dist in field (frac)')
         ylabel('Spike phase (deg)')
           title('Distance in field phase precession')
         %}
hold on

figure
title('relative to prior x')

[Nx,cx]=histcounts(x,nBinsX);
pxy=N/sum(N(:));
px=Nx/sum(Nx);
px=smooth(px);


for xi=1:nBinsX
    pYGivenX(xi,:)=pxy(xi,:)/px(xi);
end

pYGivenX=nanGaussSmooth(pYGivenX);

omarPcolor(c{1},c{2},pYGivenX',fH)
          cb=colorbar
         ylabel(cb,'Joint probability')
         colormap(gca,jet)
         title(sprintf('circ R= %.2f, p=%.3f',rho,p))
     


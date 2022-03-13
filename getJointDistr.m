function [pxySmooth,pxy]=getJointDistr(x,y,xEdges,yEdges,fH,withSmooth,cmapName)

    addPaddingToIncludeBounds=1;
    
    if(addPaddingToIncludeBounds)
        xPadAmt=range(xEdges)*0.01;
        yPadAmt=range(yEdges)*0.01;
        
        xEdges(1)=xEdges(1)-xPadAmt;
        xEdges(end)=xEdges(end)+xPadAmt;
        
        yEdges(1)=yEdges(1)-yPadAmt;
        yEdges(end)=yEdges(end)+yPadAmt;
    end
    
    if(~exist('cmapName'))
        cmapName=jet;
    end
    if(~exist('withSmooth'))
        withSmooth=0;
    end
    if(exist('fH','var'))
        showPlots=1;
    else
        showPlots=0;
    end

    %assumes y is circ variable
    nBinsX=length(xEdges)-1;
    nBinsY=length(yEdges)-1;

    %[N,c]=hist3([x(:) y(:)],[nBinsX nBinsY]);
    [N,c]=hist3([x(:) y(:)],'Edges',{xEdges,yEdges});
    N=N((1:(end-1)),1:(end-1)); %extra bin, see hist3 documentation
    c{1}=c{1}(1:(end-1));
    c{2}=c{2}(1:(end-1));
    
    pxy=N/sum(N(:));
    
    if(withSmooth)
        %pxySmooth=nanGaussSmoothPhaseProb(pxy);
        pxySmooth=nanGaussSmoothPhaseProbExtra(pxy);
    else
         pxySmooth=pxy;
    end
    pxySmooth=pxySmooth/sum(pxySmooth(:));
    
    if(showPlots)
        if(exist('fH','var'))
         omarPcolor(c{1},c{2},pxySmooth',fH)
        else
           omarPcolor(c{1},c{2},pxySmooth') 
        end
              cb=colorbar;
             ylabel(cb,'Joint probability')
             colormap(gca,cmapName)
         
    end
         %title(sprintf('circ R= %.2f, p=%.3f',rho,p))
         %{
         xlabel('Dist in field (frac)')
         ylabel('Spike phase (deg)')
           title('Distance in field phase precession')
         %}
hold on

%{
fHGivenPos=figure
title('relative to prior x')

[Nx,cx]=histcounts(x,xEdges);

px=Nx/sum(Nx);
px=smooth(px);


for xi=1:nBinsX
    pYGivenX(xi,:)=pxy(xi,:)/px(xi);
end

pYGivenX=nanGaussSmooth(pYGivenX);

        omarPcolor(c{1},c{2},pYGivenX',fHGivenPos)
  
          cb=colorbar
         ylabel(cb,'Joint probability')
         colormap(gca,jet)
         %title(sprintf('circ R= %.2f, p=%.3f',rho,p))
%}
     


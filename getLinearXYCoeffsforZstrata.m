function [] = getLinearXYCoeffsforZstrata(x,y,z)
%divides matrix into 10 z strata and computes linear relation between
%x and y for each

   x(x>1)=NaN;
        y(y>1)=NaN;
        z(z>1)=NaN;
        anyNaNrow=(isnan(x) | isnan(y) | isnan(z));
        x(anyNaNrow)=[];
        y(anyNaNrow)=[];
        z(anyNaNrow)=[];
        
numZStrata=26;

zStrataEdges=linspace(0,1,numZStrata);
zStrataMids=edgesToBins(zStrataEdges);
zStrataWidth=median(diff(zStrataMids));

figure

rPerStrata=NaN(size(zStrataMids));
for zsi=1:length(zStrataMids)
    currZstart=zStrataMids(zsi)-zStrataWidth/2;
    currZend=zStrataMids(zsi)+zStrataWidth/2;

    inStrataIdxes=z>=currZstart & z<currZend;

    currStrataX=x(inStrataIdxes);
    currStrataY=y(inStrataIdxes);
    subplot(5,5,zsi)
    plot(currStrataX,currStrataY,'k.')
    xlabel('Time in field (field frac)')
    ylabel('Spike phase (cycle frac)')
    
    %{
    cMat=corrcoef([currStrataX(:) currStrataY(:)]);
    corrCoeff=cMat(1,2);
    %}
    [ rho,p,s,b ] = kempter_lincirc( currStrataX,currStrataY*2*pi);
    corrCoeff=rho;
    
    title({sprintf('Dist from to %.2f to %.2f',abs(currZstart),currZend),sprintf('circlin R=%.3f',corrCoeff)})
    rPerStrata(zsi)=corrCoeff;
    daspect([1 1 1])
end
uberTitle(sprintf('Average phase vs time, controlling for distance in field: circ R=%.2f+/-%.2f',nanmean(rPerStrata),nanstd(rPerStrata)))
maxFig
setFigFontTo(13)
saveas(gcf,'timeVsPhaseLinearStats.tif')
disp('')




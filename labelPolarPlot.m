function labelPolarPlot(thetaLabel,rhoLabel,rmax)

    %rmax = 2;
    %hax = polaraxes('RLim', [0 rmax]);
    text(0, rmax/2, rhoLabel, 'horiz', 'center', 'vert', 'top', 'rotation', 0);
    text(pi/4, rmax*1.2, thetaLabel, 'horiz', 'center', 'rotation', -45);


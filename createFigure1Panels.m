%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description: createFigure1Panels
%
%
%Preconditions:
%
%
%Effects:
%
%
%Author: Tibin John, tibintj@umich.edu
%Project directory name: /nfs/turbo/lsa-ojahmed/tibin/spikeDynamicsAnalysisTibin 
%Created on 2018-08-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=30000;

plotPanelsSeparately=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task #1
%Description:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')
cellPropDir='/nfs/turbo/lsa-ojahmed/tibin/tibinCellSzClassificationFiles'

figure
x0=0.8;
y0=0;
width=0.2;
height=0.2;
%[axCa] = subplotCartesian(x0,y0,width,height);
pos=[x0 (y0-height) width height];
[axCa]=subplot('Position',pos)
%subplot(2,1,1)
if(plotPanelsSeparately)
	figure
	plotAllWaveformsInClass(0)
	ylim([-1 1])
	title('RS')
	print('-r600',gcf,'ClassifiedRS_AvgWaveforms','-depsc','-tiff')
else
	plotAllWaveformsInClass(0)
end


%subplot(2,1,2)
title('RS')

x0=0.8;
y0=0.325;
width=0.2;
height=0.2;
%[axCb] = subplotCartesian(x0,y0,width,height);
pos=[x0 (y0-height) width height];
[axCb]=subplot('Position',pos)

if(plotPanelsSeparately)
	figure
	plotAllWaveformsInClass(1)
	title('FS')
	ylim([-1 1])
	print('-r600',gcf,'ClassifiedFS_AvgWaveforms','-depsc','-tiff')
else
	plotAllWaveformsInClass(1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task #2
%Description:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')

x0=0.8;
y0=0.75;
width=0.18;
height=0.18;
pos=[x0 (y0-height) width height];
[axE]=subplot('Position',pos)
%[axE] = subplotCartesian(x0,y0,width,height);

x0=0;
y0=0;
width=0.45;
height=width*4/3;
%height=0.65;
%[axD] = subplotCartesian(x0,y0,width,height);
pos=[x0 (y0-height) width height];
[axD]=subplot('Position',pos)

if(plotPanelsSeparately)
	classifyFromCellProp_Tibin_Plot(cellPropDir);
else
	classifyFromCellProp_Tibin_Plot(cellPropDir,axE,axD);
end

%{
if(plotPanelsSeparately)
	figD=figure;
	copyobj(axD,figD);
	axD.OuterPosition=[ 0 0 1 1 ];
	print('-r600',gcf,'ClassifiedUnit_Mdistances','-depsc','-tiff')

	close(figD)
	
	figE=figure;
	copyobj(axE,figE);

	axE.OuterPosition=[ 0 0 1 1 ];
	print('-r600',gcf,'ClassificationWaveformFeatureSpace','-depsc','-tiff')
	close(figE)
end
%}

%axes(axE)
%ylim([0 75])
%xlim([0 75])

%axes(axD)
%xlim([0.3 0.7]) %TO DO: identify and take out noise cell! done in classifyFrom... file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task #3
%Description:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing output.........')

setFigFontTo(15)
maxFig

axes(axE)
setFontTo(8)

if(~plotPanelsSeparately)
	saveas(gcf,'humanUnitClassificationFigure.tif')
end

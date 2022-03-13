function [] =maxFigManual(sizeFact)
h = gcf; % Current figure handle
if(~exist('sizeFact'))
	sizeFact=6;
end
xLength=18;
%yLength=6;
yLength=9;

set(h,'Resize','off');
set(h,'PaperPositionMode','manual');
set(h,'PaperPosition',[0 0 xLength yLength]*sizeFact);
set(h,'PaperUnits','centimeters');
set(h,'PaperSize',[xLength yLength]*sizeFact); % IEEE columnwidth = 9cm
set(h,'Position',[0 0 xLength yLength]*sizeFact);

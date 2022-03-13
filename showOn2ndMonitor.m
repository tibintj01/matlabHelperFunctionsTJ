    MP=get(0,'MonitorPositions');
if size(MP,1)>1
    pos=get(gcf,'Position');
    pause(0.01); % this seems sometimes necessary on a Mac
    set(gcf,'Position',[pos(1,2)+MP(2,1:2) pos(3:4)]);
end
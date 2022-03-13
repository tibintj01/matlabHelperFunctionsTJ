function [] = plot3Dvectors(pts,origin,colorStr)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if(size(pts,2)<3)
    pts=[pts zeros(size(pts,1),3-size(pts,2))];
end
    
if(~exist('origin','var'))
    origin = [0,0,0];
end

numPts=size(pts,1);

if(~exist('colorStr','var'))
    colorStr=repelem('k',numPts);
end

for i=1:numPts
    plot3(origin(1),origin(2),origin(3),'r*')
    hold on
    plot3([origin(1) pts(i,1)],[origin(2) pts(i,2)],[origin(3) pts(i,3)],'LineWidth',5,'Color',colorStr(i));
   
end
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
%set(gca,'CameraPosition',[2 2 2]);
end


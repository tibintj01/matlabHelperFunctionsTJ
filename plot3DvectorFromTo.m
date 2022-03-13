function [] = plot3DvectorFromTo(origins,destinations,colorStr)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if(nargin==1)
    destinations=origins;
    origins = repmat([0,0,0],size(destinations,1),1);
end

if(~exist('origins','var'))
    origins = repmat([0,0,0],size(destinations,1),1);
end

numPts=size(destinations,1);

if(~exist('colorStr','var'))
    colorStr=repelem('k',numPts);
end

for i=1:numPts
    plot3(origins(i,1),origins(i,2),origins(i,3),'r*')
    hold on
    plot3([origins(i,1) destinations(i,1)],[origins(i,2) destinations(i,2)],[origins(i,3) destinations(i,3)],'LineWidth',5,'Color',colorStr(i));
   
end
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
%set(gca,'CameraPosition',[2 2 2]);
end


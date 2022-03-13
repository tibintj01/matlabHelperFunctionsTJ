function [degreeOnCone]=getConeDegPredForXandTvals(xVal,tVal)
%return deg (0, 360) for x (0,1) and t (0,1)

showPlots=1;
showPlots=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create cone with vertex at origin and increasing radius to 1 at z=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;
maxZ=0.25;
maxZ=0.45;
maxZ=0.35;
maxZ=0.5;
maxZ=0.45;
maxZ=0.35;
maxZ=0.33;
maxZ=0.335;


%maxZ=0.28;

%r = linspace(0, 1, N);
r = linspace(0, maxZ, N);
[X,Y,Z] = cylinder(r, N); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert to vector format (each point is a 3D column vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyzPoints=[X(:)'; Y(:)'; Z(:)']; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rotate 90 degrees towards x axis (around y axis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotTowardsX=roty(90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rotate 90 degrees towards y axis (around z axis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotTowardsY=rotx(-45);


compositeRotation=rotTowardsX*rotTowardsY;
%compositeRotation=rotTowardsY;

xyzPoints=compositeRotation*xyzPoints;

X=reshape(xyzPoints(1,:),size(X));
Y=reshape(xyzPoints(2,:),size(Y));
Z=reshape(xyzPoints(3,:),size(Z));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalize Z so that cone goes through point (1,1,1) [normalized]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Z=Z/2;
spacetimeStrechFact=1.5;
X=X*spacetimeStrechFact;
Y=Y*spacetimeStrechFact;


%{
[azimuth,elevation,r] = cart2sph(X,Y,Z);
azimuth=azimuth+pi/6;

%elevation=elevation+pi/4;

azimuth=mod(azimuth,2*pi);
elevation=mod(elevation,2*pi);
%elevation=zeros(size(elevation));

[X,Y,Z] = sph2cart(azimuth,elevation,r);


p=0;
 X = X;
 Y = Y*cos(p)-Z*sin(p);
 Z = Y*sin(p)+Z*cos(p);
%}
%%
Z(Z>0)=NaN;
Z=scaledata(Z,0,360);

tolLev=0.01;

closeXIdxes=abs(X-xVal)<tolLev;
closeTIdxes=abs(Y-tVal)<tolLev;

%{
figure;
subplot(1,3,1); imagesc(closeXIdxes);daspect([1 1 1])
subplot(1,3,2); imagesc(closeTIdxes);daspect([1 1 1])
subplot(1,3,3); imagesc(closeXIdxes & closeTIdxes & ~isnan(Z));daspect([1 1 1])
%}


correspondingZIdx=find(closeXIdxes & closeTIdxes & ~isnan(Z));

if(isempty(correspondingZIdx))
    degreeOnCone=NaN;
else
    correspondingZIdx=correspondingZIdx(1);

    degreeOnCone=Z(correspondingZIdx);
end


%{
XInt=round(X*N*10);
YInt=round(Y*N*10);

xValInt=round(xVal*N*10);
tValInt=round(tVal*N*10);
%}
disp('')


%{
for yi=1:size(X,1)
    for xi=1:size(X,2)
        
    end
end

[~,linIDx]=min(abs(X(:)-xVal))
Y-tVal

[~, index] = ismember( [xVal, tVal], [X(:) Y(:)], 'rows' );
%}

%{
xIdx=round(xVal*N);
tIdx=round(tVal*N);
xVal=X(xIdx,tIdx);

%degreeOnCone=Z()

%}
if(showPlots)
    figure
    subplot(1,2,1)
    %h = surf(X, Y, Z);
    h = surf(X, Y, Z);
    %rotate(h, [-1 1 0], 90);
    daspect([1 1 1])
    %xlim([0 1]); ylim([0 1])
    %caxis([0 maxZ])
    %zlim([0 maxZ])
    subplot(1,2,2)
    figure
    %coneMax=max(Z(:));
    %Z=scaledata(Z,0,1);

    %{
    Z=1-Z;

    Z=Z-min(Z(:));
    Z=Z/0.3;

    %}


    subplot(1,2,1)
    %[C,h]=contour(X, Y, Z,100);
    [C,h]=contour(X, Y, Z,50);
    h.LineWidth=3;

    cb=colorbar
    colormap(jet)




    caxis([0 360])

    xlim([0 1]); ylim([0 1])
    daspect([1 1 1])

    hold on
    plot([0 0.7],[0 1],'k--','LineWidth',3)
    plot([0 1],[0 0.7],'k--','LineWidth',3)


    subplot(1,2,2)

    surf(X, Y, Z);

    colorbar
    colormap(jet)

    xlim([0 1]); ylim([0 1])
    daspect([1 1 360])

    view( -35.5000,63.6000)
end






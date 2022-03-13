close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
coneMax=max(Z(:));
%Z=scaledata(Z,0,1);

%{
Z=1-Z;

Z=Z-min(Z(:));
Z=Z/0.3;

%}
Z(Z>0)=NaN;
Z=scaledata(Z,0,360);

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






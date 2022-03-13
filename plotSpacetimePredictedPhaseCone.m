function [Zq] = plotSpacetimePredictedPhaseCone(xyPairs,surfPlotStr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create cone with vertex at origin and increasing radius to 1 at z=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
showPlots=1;
if(exist('xyPairs','var'))
    showPlots=0;
end

dispVersion=1;

N = 1000;

if(exist('surfPlotStr','var'))
    N=100;
    maxZ=1;
end

%{
maxZ=0.25;
maxZ=0.45;
maxZ=0.35;
maxZ=0.5;
maxZ=0.45;
maxZ=0.35;
maxZ=0.33;
%maxZ=0.28;
%}

maxZ=0.335;
maxZ=0.28;
maxZ=0.33;
%maxZ=0.335;




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
if(~dispVersion)
rotTowardsY=rotx(-51);
end
%rotTowardsY=rotx(-53);


compositeRotation=rotTowardsX*rotTowardsY;
%compositeRotation=rotTowardsY;

xyzPoints=compositeRotation*xyzPoints;

originalXsize=size(X);
originalYsize=size(Y);
originalZsize=size(Z);

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


%X=X+0.02; %shift to right
if(~dispVersion)
X=X-0.075*1.1; %shift to left
Y=Y-0.075*1.1; %shift down
end

%goodCoordinates=X<=1 & Y<=1 & Z>=0;
%X(X>1)=NaN;
%Y(Y>1)=NaN;

Z(Z<0)=NaN;


Z=1-Z/maxZ;

if(exist('surfPlotStr','var'))
    %figure(fH)
    surf(X, Y, Z)
    colormap(gca,'jet')
    daspect([1 1 1])
    caxis([0 1])
    zlim([0 1])
    box off
    grid off
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get closest z value for each xy pair in input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(showPlots==0)
    numPts=size(xyPairs,1);
    Zq=NaN(numPts,1);
    tol=7.0605e-04;
    for i=1:numPts
        currXq=xyPairs(i,1);
        currYq=xyPairs(i,2);
        %Zq=interp2(X,Y,Z,currXq,currYq);

        %[closestX,closestXid]=min(abs(X(:,1)-currXq));
        %[closestY,closestYid]=min(abs(Y(:,1)-currYq));
        %closestXid=round(currXq*N/spacetimeStrechFact);
        %closestYid=round(currYq*N/spacetimeStrechFact);

        validIDs=find(abs(X-currXq)<tol & abs(Y-currYq)<tol & ~isnan(Z));
        if(~isempty(validIDs))
            Zq(i)=Z(validIDs(1));
        end
        %Zq(i)=Z(closestYid,closestXid);

        %{
        zInd=sub2ind(size(Z),closestYid,closestXid);

        linearizedZ=Z(:);
        currZq=linearizedZ(zInd);
        Zq(i)=currZq;
        %}



    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot contours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(showPlots)
    [C,h]=contour(X, Y, Z,50);
    colormap(gca,jet)
    caxis([0 1])
    h.LineWidth=3;
    cb=colorbar;
    ylabel(cb,'Predicted theta phase (cycle frac)')

    xlim([0 1]); ylim([0 1])
    daspect([1 1 1])
           xlabel('Distance in field (frac)')
   ylabel('Time in field (frac)')
   title('Spacetime cone model')
   box off
   
end

%{
hold on
plot([0 0.7],[0 1],'k--','LineWidth',3)
plot([0 1],[0 0.7],'k--','LineWidth',3)
%}






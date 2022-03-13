function isAboveEllipse=isAboveOrOnEllipse(pts,ellipseParams,fH)
%returns whether each 3D pt in pts is above ellipse of given parameters
    numPts=size(pts,1);
    
    isAboveEllipse=zeros(numPts,1);
    
    C=ellipseParams.centerPt;
    
    a=ellipseParams.majorAxisLength;
    b=ellipseParams.minorAxisLength;
    
    U=ellipseParams.majorAxisUnitVector;
    V=ellipseParams.minorAxisUnitVector;
    
    angles=0:0.00005:2*pi;
    
    for ai=1:length(angles);
        currAngle=angles(ai);
        xE(ai)=C(1)+a*cos(currAngle)*U(1)+b*sin(currAngle)*V(1);
        yE(ai)=C(2)+a*cos(currAngle)*U(2)+b*sin(currAngle)*V(2);
        zE(ai)=C(3)+a*cos(currAngle)*U(3)+b*sin(currAngle)*V(3);
    end
    
    if(exist('fH','var'))
        figure(fH); plot3(xE(1:1000:end),yE(1:1000:end),zE(1:1000:end),'Color',ellipseParams.colorStr)
        hold on
    end
    
    for i=1:numPts
        currPt=pts(i,:);
        
        %get angle corresponding to current point
        currAngle=atan2(currPt(2),currPt(1));
        [~,currAi]=min(abs(angles-currAngle));
        if(zE(currAi)<=currPt(3))
            isAboveEllipse(i)=1;
        end
        
    end
    disp('here')
    
end


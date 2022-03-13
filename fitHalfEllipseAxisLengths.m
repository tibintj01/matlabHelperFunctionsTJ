function  [bestW,bestH]=fitHalfEllipseAxisLengths(ellipseCenter,x,y,fixedA2)
%finds best semi minor and major axis lengths to fit data that resembles
%half-ellipse with given center
    useFixed=1;
    if(~exist('fixedA2'))
        fixedA2=1;
        useFixed=0;
    end

    numPts=length(x);
    
    X0=ellipseCenter(1);
    Y0=ellipseCenter(2);
    
    %t = linspace(0,pi,numPts);
    t=atan2(y-Y0,x-X0);
    t=t(:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model:     
    %a=[0 0];
    %x=a(1) * cos(t);
    %y=a(2) * sin(t);
    %find least squares a: (W*a-x)
    %W is (2*numPts( x 2 matrix
    %a is 2 x 1
    %x is 2*numPts x 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W=[cos(t(:)) zeros(numPts,1) ; ...
        zeros(numPts,1) fixedA2*sin(t(:))];
    
    af= W \ [(x-X0);(y-Y0)];
    
      
    if(useFixed==1)
        af(2)=fixedA2;
    end
    xnew = X0+af(1)*cos(t);
    ynew = Y0+af(2)*sin(t);
  
    bestW=af(1);
    bestH=af(2);
   

    figure
    plot(x,y,'.')
    hold on
    plot(xnew,ynew,'ro','LineWidth',3)

end


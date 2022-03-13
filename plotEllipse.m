function [ellipseX,ellipseY] = plotEllipse(ellipseStruct)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
showPlots=0;
numPts=1000;
    phi=ellipseStruct.phi;
    cos_phi=cos(phi);
    sin_phi=sin(phi);
    
    a=ellipseStruct.a;
    b=ellipseStruct.b;
    X0=ellipseStruct.X0;
    Y0=ellipseStruct.Y0;

% rotation matrix to rotate the axes with respect to an angle phi
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    
    
    
    % the axes
    ver_line        = [ [X0 X0]; Y0+b*[-1 1] ];
    horz_line       = [ X0+a*[-1 1]; [Y0 Y0] ];
    new_ver_line    = R*ver_line;
    new_horz_line   = R*horz_line;
    
    % the ellipse
    theta_r         = linspace(0,2*pi,numPts);
    ellipse_x_r     = X0 + a*cos( theta_r );
    ellipse_y_r     = Y0 + b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    
    % draw
    %hold_state = get( axis_handle,'NextPlot' );
    %set( axis_handle,'NextPlot','add' );
    if(showPlots)
        plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
        plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
        plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r','LineWidth',5);
    end
    
    ellipseX=rotated_ellipse(1,:);
    ellipseY=rotated_ellipse(2,:);
end


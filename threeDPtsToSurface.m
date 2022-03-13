function [zs,avgs,counts] = threeDPtsToSurface(data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [ux, ~, xidx] = unique(data(1,:));
    [uy, ~, yidx] = unique(data(2,:));

    %count the number of points at each unique x/y combination
    counts = accumarray([xidx(:), yidx(:)], 1);  
    %average the z that fall into each unique x/y combination
    avgs = accumarray([xidx(:), yidx(:)], data(3,:).');
    %create a list of the z that fall into each unique x/y combination
    %zs = accumarray([xidx(:), yidx(:)], data(3,:).', [], @(V) {V}, {});
    zs = accumarray([xidx(:), yidx(:)], data(3,:).', [], @(V) {V}, {NaN});
end


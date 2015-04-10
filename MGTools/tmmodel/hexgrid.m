function [X, Y, dist] = hexgrid(npoints)
% HEXGRID generate coordinates of grids in range [-1,1] x [-1,1] with
% hexagon shape
%
% USAGE  [X,Y] = HEXGRID(NPOINTS) would generate a grid around NPOINTS and
%        distribute as uniform as possible
%
% MooGu Z. <hzhu@case.edu>
% Apr 9, 2015

% estimate number of rows and columns
rt = sqrt(2*npoints-1);
% residual of estimation
r = rt - floor(rt);
% completely square number
if r == 0 
    nrow = rt;
    ncol = rt;
% otherwise, need arrange the points to find a combination of nrow and ncol
else
    th = (rt + 1/2) - sqrt(rt^2 + 1/4);
    % make nrow less than ncol
    if r <= th
        nrow = floor(rt);
        ncol = ceil(rt);
    else
        nrow = ceil(rt);
        ncol = ceil(rt);
    end
end

% calculate number of point under current arrangement
npts = nrow * floor(ncol/2) + rem(ncol,2) * ceil(nrow/2);
% interval between rows
rint = 2 / (nrow - 1);
% interval between columns
cint = 2 / (ncol - 1);

% generate first half of points
[Xa, Ya] = meshgrid( ...
    linspace(-1, cond(rem(ncol,2), 1, 1-cint), ceil(ncol/2)), ...
    linspace(1, cond(rem(nrow,2), -1, rint-1), ceil(nrow/2)));
% generate second half of points
[Xb, Yb] = meshgrid( ...
    linspace(cint-1, cond(rem(ncol,2), 1-cint, 1), floor(ncol/2)), ...
    linspace(1-rint, cond(rem(nrow,2), rint-1, -1), floor(nrow/2)));
% combine points
X = [Xa(:);Xb(:)];
Y = [Ya(:);Yb(:)];

% remove periphery points to fit requirement of NPOINTS, if necessary
if npts > npoints
    [~,index] = sort(X.^2 + Y.^2, 'ascend');
    % preserve NPOINTS points close to [0,0]
    X = X(index(1:npoints));
    Y = Y(index(1:npoints));
end

% calculate minimum distance between each points
dist = min([2 * rint, 2 * cint, sqrt(rint^2 + cint^2)]);

end

function aorb = cond(bool, a, b)
% this helper function work as the operator 'bool ? a : b' in C, here you
% use cond(bool, a, b)

if bool
    aorb = a;
else
    aorb = b;
end

end

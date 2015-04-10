function T = tbaseGen(orient,freq,res)
% TBASEGEN generate transform base function (set)
%
% USAGE T = tbaseGen(ORIENT,FREQ=3,RES=256) generate
%       transformation base with ORIENT (angle of
%       phase gradient), frequency FREQ in RESxRES
%       image. Where, both ORIENT and FREQ can be
%       to generate a set of functions. ORIENT is in
%       unit degree.
%
% MooGu Z. <hzhu@case.edu>
% Apr 9, 2015

% set default values
if ~exist('freq','var'), freq = 3; end
if ~exist('res','var'), res = 256; end

% number of base functions
nbase = numel(orient) * numel(freq);
% initialize base functions
T = zeros(res,res,nbase);
% calculate coordinates
[X,Y] = meshgrid(linspace(-pi,pi,res),linspace(pi,-pi,res));
% generate bases one by one
ibase = 1;
for i = 1 : numel(orient)
    theta = pi * orient(i) / 180;
    for j = 1 : numel(freq)
        f = freq(j);
        T(:, :, ibase) = ...
            wrapToPi(f * (cos(theta) * X + sin(theta) * Y));
        ibase = ibase + 1;
    end
end

end

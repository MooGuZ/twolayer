function bases = tbaseGen(type, paramA,paramB,res)
% TBASEGEN generate transform base function (set)
%
% BASES = tbaseGen('translate', ORIENT, FREQ, RES=256) generate
%   transformation for translation with ORIENT in degree (from 0 to 360)
%   and FREQ, which corresponding to the number of cycles in the range of
%   current base. Both ORIENT and FREQ can be a list. If this is the case,
%   all combination between them would produce a corresponding base in the
%   result. RES define the resolution of base function as a discrete image,
%   and default to 256x256.
%
% BASES = tbaseGen('scale', CENTER, FREQ, RES=256) generate transformation
%   bases for scaling with center at CENTER (range from [-1,-1] to [1,1])
%   and frequency FREQ, which represent the number of cycle in a distance
%   unit (here is the half width of the base). Both CENTER and FREQ can be
%   a list. Then all combination between them would produce one base in the
%   result. RES define the resolution of base function as a discrete image,
%   and default to 256x256.
%
% BASES = tbaseGen('rotate', CENTER, NCYCLE, RES=256) generate
%   transformation bases for rotating, which is centered at CENTER (range
%   from [-1,-1] to [1,1]) and number of cycles, NCYCLE, in a round. Both
%   CENTER and NCYCLE can be a list. Then all combination between them 
%   would produce one base in the result. RES define the resolution of base 
%   function as a discrete image, and default to 256x256.
%
%   see also mbaseGen.
%
% MooGu Z. <hzhu@case.edu>
% Apr 9, 2015

% set default values
if ~exist('res', 'var'), res = 256; end
if numel(res) == 1, res = [res, res]; end

switch lower(type)
    case {'translate', 'translation', 'translating', 'shift', 'shifting'}
        bases = translateBases(paramA, paramB, res);
        
    case {'rotation', 'rotate', 'rotating'}
        bases = rotateBases(paramA, paramB, res);
        
    case {'scale', 'scaling'}
        bases = scaleBases(paramA, paramB, res);
        
    otherwise
        error('Unrecognized transformation type : %s\n', type)
end

end

function bases = translateBases(orient, freq, res)
% this helper function generate transformation bases for translation.

% number of base functions
nbase = numel(orient) * numel(freq);
% initialize base functions
bases = zeros([res,nbase]);
% calculate coordinates
[X,Y] = meshgrid(linspace(-pi,pi,res(2)),linspace(pi,-pi,res(1)));
% generate bases one by one
for i = 1 : numel(orient)
    theta = pi * orient(i) / 180;
    for j = 1 : numel(freq)
        bases(:, :, (i-1) * numel(freq) + j) = ...
            wrapToPi(freq(j) * (cos(theta) * X + sin(theta) * Y));
    end
end

end

function bases = rotateBases(center, freq, res)
% this helper function generate transformationa bases for rotation.

% number of base functions needed to generate
nbase = size(center,1) * numel(freq);
% initialize base funcitons
bases = zeros([res,nbase]);
% create coordinates
[X,Y] = meshgrid(linspace(-1,1,res(2)),linspace(1,-1,res(1)));
% generate bases
for i = 1 : size(center, 1)
    for j = 1 : numel(freq)
        bases(:, :, (i-1) * numel(freq) + j) = wrapToPi( ...
            angle((X - center(i,1)) + 1j*(Y - center(i,2))) * freq(j));
    end
end

end

function bases = scaleBases(center, ncycle, res)
% this helper function generate tranaformation bases for scaling.

% number of base functions needed to generate
nbase = size(center,1) * numel(ncycle);
% initialize base funcitons
bases = zeros([res,nbase]);
% create coordinates
[X,Y] = meshgrid(linspace(-pi,pi,res(2)),linspace(pi,-pi,res(1)));
% generate bases
for i = 1 : size(center, 1)
    for j = 1 : numel(ncycle)
        bases(:, :, (i-1) * numel(ncycle) + j) = ...
            reshape(wrapToPi(sqrt(sum( ...
                bsxfun(@minus,[X(:),Y(:)],center(i,:)) .^ 2, ...
            2)) * ncycle(j)), res, res);
    end
end

end

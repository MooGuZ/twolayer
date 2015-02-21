function T = tbaseGen(ovec,freq,res)
% TBASEGEN would generate a simple transform base function

% check input argument
assert(numel(ovec) == 2, ...
    'orientation vector should have and only have two elements.')

% set default values
if ~exist('freq','var'), freq = 3; end
if ~exist('res','var'), res = 256; end

% calculate coordinates
[X,Y] = meshgrid(linspace(-freq*pi,freq*pi,res), ...
    linspace(freq*pi,-freq*pi,res));
% normalize orientation vector
ovec = ovec / norm(ovec);
% generate transform base on specific orientation
T = wrapToPi(ovec(1) * X + ovec(2) * Y);

end

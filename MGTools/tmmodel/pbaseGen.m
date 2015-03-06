function P = pbaseGen(cvec,radius,res)
% PBASEGEN would generate a Gaussian mask as a simple pattern base function

% check input arguments
assert(numel(cvec) == 2, ...
    'center vector should have and only have two elements.');

% set default values
if ~exist('radius','var'), radius = 4; end
if ~exist('res','var'), res = 8 * radius; end

% calculate coordinates
[X,Y] = meshgrid(linspace(-1,1,res),linspace(1,-1,res));
% calculate sigma for Gaussian function
sigma = radius / res;
% generate pattern base
P = exp(-((X - cvec(1)).^2 + (Y - cvec(2)).^2) / (2*sigma^2));
P = P / max(P(:));

end
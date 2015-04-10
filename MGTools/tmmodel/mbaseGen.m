function M = mbaseGen(cvec,radius,res)
% PBASEGEN would generate a Gaussian mask as a simple pattern base function

% check input arguments
assert(size(cvec,1) == 2, ...
    'center vector should have and only have two elements.');

% set default values
if ~exist('radius','var'), radius = 4; end
if ~exist('res','var'), res = 8 * radius(1); end

% calculate number of base function
nbase = size(cvec,2) * numel(radius);

% initialize bases
M = zeros(res,res,nbase);
% calculate coordinates
[X,Y] = meshgrid(linspace(-1,1,res),linspace(1,-1,res));
% generate base function one by one
for i = 1 : numel(radius)
    % calculate sigma for Gaussian function
    sigma = radius(i) / res;
    for j = 1 : size(cvec,2)
        % generate mask base
        buffer = exp(-((X - cvec(1,j)).^2 + (Y - cvec(2,j)).^2) ...
            / (2*sigma^2));
        M(:,:,j) = buffer / max(buffer(:));
    end
end

end
function anim = animGen(shape,trans,varargin)
% ANIMGEN is a generator for several kind of animation.
%
% MooGu Z. <hzhu@case.edu>
% Version 0.1 - Jan 20, 2015
% Version 0.5 - Jan 29, 2015

% define available shape and transformation set
shapeSet = {'circle','edge','sine','triangle','square','polygon'};
transSet = {'move','rotate','scale','movescale','moverotate', ...
    'rotatescale','mix'};
% special transformation sets
moveTransSet   = {'move','moverotate','movescale','mix'};
rotateTransSet = {'rotate','moverotate','rotatescale','mix'};
scaleTransSet  = {'scale','movescale','rotatescale','mix'};
% initialize input parser
p = inputParser;
% set input parameters
p.addRequired('shape',@(x) any(validatestring(lower(x),shapeSet)));
p.addRequired('trans',@(x) any(validatestring(lower(x),transSet)));
p.addParamValue('frmsize',[32,32], ...
    @(x) isnumeric(x) && isreal(x) && all((x-floor(x))==0) ...
        && numel(x) <= 2);
p.addParamValue('frmnumber',24, ...
    @(x) isscalar(x) && isnumeric(x) && isreal(x));
p.addParamValue('objsize',10, ...
    @(x) isnumeric(x) && isreal(x) && x > 0);
p.addParamValue('objorient',0, ...
    @(x) isscalar(x) && isnumeric(x) && isreal(x));
p.addParamValue('objposition',nan, ...
    @(x) any(isnan(x)) || (isnumeric(x) && isreal(x) && numel(x) == 2));
p.addParamValue('smoothness',1, ...
    @(x) isscalar(x) && isnumeric(x) && isreal(x) && x >= 0);
p.addParamValue('mvdirection',0, ...
    @(x) isscalar(x) && isnumeric(x) && isreal(x));
p.addParamValue('edgenumber',5, ...
    @(x) isscalar(x) && isnumeric(x) && isreal(x) && x - floor(x) == 0);
% parsering input arguments
p.parse(shape,trans,varargin{:});
% interpret parsering results
frmsz     = p.Results.frmsize;
nfrm      = p.Results.frmnumber;
objsz     = p.Results.objsize;
objorient = p.Results.objorient;
objpos    = p.Results.objposition;
tzwide    = p.Results.smoothness * (2 / objsz);
mvdirct   = p.Results.mvdirection;
nedge     = p.Results.edgenumber;
% initialize transform flags
moving   = any(strcmpi(trans,moveTransSet));
rotating = any(strcmpi(trans,rotateTransSet));
scaling  = any(strcmpi(trans,scaleTransSet));

% parameter standardization
if isscalar(frmsz), frmsz = [frmsz,frmsz]; end
nfrm = round(nfrm);

% initialize coordinates
[X,Y]  = meshgrid(linspace(-1,1,frmsz(1)),linspace(1,-1,frmsz(2)));
xscale = frmsz(1) / objsz;
yscale = frmsz(2) / objsz;
Z      = complex(xscale * X, yscale * Y);
% set start position for the object
if ~isnan(objpos), Z = Z - objpos * [xscale;1j*yscale]; end
% adjust coordinates to start status according to transformation type
if moving
    % move vector (in complex number form)
    mV = exp(1j*mvdirct);
    % shift to farthest start point, if have not specified object point
    Z  = Z - min(real(Z(:) * conj(mV))) * mV;
    % rotate move vector to compansate the effects of rotation coordinate
    mV = mV * exp(-1j*objorient);
end
% rotate coordinates to fit objec orientation
Z  = Z * exp(-1j*objorient);

% initialize animation data
anim = zeros(prod(frmsz),nfrm);
% calculate step size for each transformation
if nfrm > 1
    if moving
        mvstep = max(real(Z(:) * conj(mV))) / nfrm;
    end
    rtstep = 2 * pi / nfrm;
    scstep = (max(abs(Z(:))))^(1/(nfrm-1));
end
% transforming frame by frame
for f = 1 : nfrm
    % generate current frame accordint to the shape
    switch lower(shape)
        case {'circle'}
            curfrm = circle(Z,@(M) bfunc(M,tzwide));
        case {'edge'}
            curfrm = edge(Z,@(M) bfunc(M,tzwide));
        case {'sine'}
            curfrm = edge(Z,@(M) bfunc(M,1));
        case {'triangle'}
            curfrm = polygon(Z,3,@(M) bfunc(M,tzwide));
        case {'square'}
            curfrm = polygon(Z,4,@(M) bfunc(M,tzwide));
        case {'polygon'}
            curfrm = polygon(Z,nedge,@(M) bfunc(M,tzwide));
        otherwise
            error('Shape does not defined!');
    end
    % record frame in animation
    anim(:,f) = curfrm(:);
    % transforming to next status
    if moving, Z = Z - mvstep * mV; end
    if rotating
        Z  = Z * exp(-1j*rtstep);
        if moving
            mV = mV * exp(-1j*rtstep);
        end
    end
    if scaling
        Z = Z / scstep;
        tzwide = tzwide / scstep;
        if moving
            mvstep = mvstep / scstep;
        end
    end
end   
            
end

function I = circle(Z,bfunc)
I = bfunc(abs(Z) - 1);
end

function I = edge(Z,bfunc)
I = bfunc(abs(real(Z)) - 1);
end

function I = polygon(Z,n,bfunc,orient,dist)
% norm vector orientation for each edge
if ~exist('orient','var')
    orient = (0 : n-1) * (2*pi / n);
end
% distance from center to edges
if ~exist('dist','var')
    dist = cos(pi/n) * ones(1,n);
end
% initialize image
I = zeros(size(Z));
% apply each edge to the image
for i = 1 : n
    normVec = complex(cos(orient(i)),sin(orient(i)));
    I = max(I,bfunc(real(Z * conj(normVec)) - dist(i)));
end

end

function I = bfunc(M,wide)
% initialize I
I = zeros(size(M));
% calculate the map of transition zone
tzone = abs(M) < wide;
% elements locate right side of tzone assign to 1
I(M >= wide) = 1;
% elements locate left side of tzone assign to 0
I(M <= -wide) = 0;
% calculate the value of tzone elements
if wide ~= 0
    I(tzone) = (sin( (pi*M(tzone)) / (2*wide) ) + 1) / 2;
%     I(tzone) = (M(tzone) / wide + 1) / 2;
end 
end

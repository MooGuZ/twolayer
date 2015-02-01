function img = comp2img(M,varargin)
% COMP2IMG convert complex matrix into RGB image, whose saturation
% modulated by norm of complex numbers and hue modulated by phase of
% complex numbers in the matrix

% Input Parsing
p = inputParser;

p.addRequired('M',@(x) ismatrix(x)&&isnumeric(x)&&~isreal(x));
p.addOptional('method','linear',...
    @(x) any(validatestring(x,{'linear','cutoff'})));
p.addParamValue('Rmax',max(abs(M(:))),...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));
p.addParamValue('Smax',1,...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));
p.addParamValue('Brightness',0.9,...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));
p.addParamValue('Resolution',256,...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));
p.addParamValue('Threshold',0.07,...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));

p.parse(M,varargin{:});

method = p.Results.method;
Rmax   = p.Results.Rmax;
Smax   = min(p.Results.Smax,1);
Bright = min(p.Results.Brightness,1);
Res    = floor(p.Results.Resolution);
Th     = p.Results.Threshold;

% Create Mapping Functions
switch lower(method)
    case 'linear'
        S = @(r) (r/Rmax)*Smax;
    case 'cutoff'
        S = @(r) double(r>Th)*Smax;
    otherwise
        error('[COMP2IMG] Undefined Modulating Method!');
end
H = @(phi) (sin(phi/2)+1)/3;
V = @(m) Bright*ones(size(m));

% Generate Image in HSV Color Space
[m,n] = size(M);
if Res^2 < 4*m*n
    img = zeros(m,n,3);
    img(:,:,1) = H(angle(M));
    img(:,:,2) = S(abs(M));
    img(:,:,3) = V(M);
else
    scale = round(Res/sqrt(m*n));
    img = zeros(scale*m,scale*n,3);
    img(:,:,1) = kron(H(angle(M)),ones(scale));
    img(:,:,2) = kron(S(abs(M)),ones(scale));
    img(:,:,3) = kron(V(M),ones(scale));
end

% Convert Image into RGB Space
img = hsv2rgb(img);

end


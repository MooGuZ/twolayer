function img = comp2img(M,varargin)
% COMP2IMG convert complex matrix into RGB image, whose brightness
% modulated by norm of complex numbers and color modulated by phase of
% complex numbers

% Input Parsing
p = inputParser;

p.addRequired('M',@(x) ismatrix(x)&&isnumeric(x)&&~isreal(x));
p.addOptional('method','linear',...
    @(x) any(validatestring(x,{'linear','cutoff'})));
p.addParamValue('Rmax',max(abs(M(:))),...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));
p.addParamValue('Ymax',0.3,...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));
p.addParamValue('Ymin',-0.9,...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));
p.addParamValue('Resolution',512,...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));
p.addParamValue('Threshold',0.07,...
    @(x) isnumeric(x)&&isscalar(x)&&isreal(x));

p.parse(M,varargin{:});

method = p.Results.method;
Rmax   = p.Results.Rmax;
Ymax   = p.Results.Ymax;
Ymin   = p.Results.Ymin;
Res    = p.Results.Resolution;
Th     = p.Results.Threshold;

% Create Mapping Functions
switch lower(method)
    case 'linear'
        Y = @(r) ((Ymax-Ymin)/Rmax)*r + Ymin;
    case 'cutoff'
        Y = @(r) (Ymax+1)*double(r>Th) - 1;
    otherwise
        error('[COMP2IMG] Undefined Modulating Method!');
end
Cb = @(phi)  0.45*sin(phi+(37/32)*pi) + 0.5;
Cr = @(phi) -0.45*sin(phi+(27/32)*pi) + 0.5;

% Generate Image in YCbCr Color Space
[m,n] = size(M);
if Res^2 < 4*m*n
    img = zeros(m,n,3);
    img(:,:,1) = Y(abs(M));
    img(:,:,2) = Cb(angle(M));
    img(:,:,3) = Cr(angle(M));
else
    scale = round(Res/sqrt(m*n));
    img = zeros(scale*m,scale*n,3);
    img(:,:,1) = kron(Y(abs(M)),ones(scale));
    img(:,:,2) = kron(Cb(angle(M)),ones(scale));
    img(:,:,3) = kron(Cr(angle(M)),ones(scale));
end

% Convert Image into RGB Space
img = ycbcr2rgb(img);

end


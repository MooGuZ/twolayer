function steerBases(Fh,Fv,varargin)
% STEERBASES show the animation of two base functions by steering them
% through weights

omSet = {'disp','gif'};

% Input Parsing
% -------------
p = inputParser;
% set input arguments
p.addRequired('HorizontalFrequencyORBaseA', ...
    @(x) isnumeric(x) && isreal(x));
p.addRequired('VerticalFrequencyORBaseB', ...
    @(x) isnumeric(x) && isreal(x));
p.addParamValue('Method','disp', ...
    @(x) any(validatestring(lower(x),omSet)));
p.addParamValue('FileName',[datestr(now),'.gif'],@isstr);
p.addParamValue('Duration',10, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && x > 0);
p.addParamValue('Resolution',256, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) ...
    && x > 0 && floor(x) == x);
% parsing
p.parse(Fh,Fv,varargin{:});
% set up parameters according to input arguments
method = p.Results.Method;
fname  = p.Results.FileName;
nframe = ceil(p.Results.Duration * 24);
res    = p.Results.Resolution;
% clear parser
clear('p');
%--------------------------------------------------------------------------

% construct initial phase map
if isscalar(Fh)
    A = repmat(linspace(-pi,pi,res)*Fh,[res,1]); 
else
    A = Fh;
end
if isscalar(Fv)
    B = repmat(linspace(-pi,pi,res)'*Fv,[1,res]);
else
    B = Fv;
end
% make angle vector
ang = linspace(0,2*pi,nframe);
% creat animation
anim = mat2img(bsxfun(@times,exp(1j*A(:)),sin(ang)) ...
    + bsxfun(@times,exp(1j*B(:)),cos(ang)));
% display or write the animation
if exist('method','var') && strcmpi(method,'gif')
    assert(exist('fname','var'),'write into GIF file need a file name.');
    anim2gif(anim,'FileName',fname);
else
    animshow(anim);
end

end
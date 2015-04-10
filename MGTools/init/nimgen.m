function I = nimgen(ntype,varargin)
% NIMGEN is short for Noise IMage GENerator. This function generate
% different types of noise images.
%
% MooGu Z. <hzhu@case.edu>
% Jan 15, 2015 - Version 0.1
% Apr 01, 2015 - Version 0.2
%   Reformed as a general noise image generator

% supported noise type list
ntypelist = {'white','pinklike','cow'};
% spectrum mask generating function list
mfunclist = {@white, @pinklike, @cow};
% supported random variable distribution
vdistlist = {'gaussian','uniform'};
% random variable generate funciton list
vfunclist = {@randn, @rand};
% binary process needed list
bproclist = {'cow'};

% initialize input parser
p = inputParser;
% define input parameters
p.addRequired('ntype', @(x) any(validatestring(x,ntypelist)));
p.addOptional('imsize', 256, @(x) isreal(x) && (numel(x) <= 2));
p.addParamValue('corr', 1, @(x) isscalar(x) && isreal(x));
p.addParamValue('cfreq', 0, @(x) isscalar(x) && isreal(x));
p.addParamValue('dist', 'uniform', @(x) any(validatestring(x,vdistlist)));
% parsing input
p.parse(ntype,varargin{:});
% initialize parameters
ntype  = p.Results.ntype;
imsize = p.Results.imsize;
corr   = p.Results.corr;
cfreq  = p.Results.cfreq;
dist   = p.Results.dist;

% generate white noise image
I = vfunclist{strcmpi(dist,vdistlist)}(imsize(end:-1:1));

% calculate spectrum of white noise
F = fftshift(fft2(I));

% apply corresponding mask to get required noise spectrum
F = F .* mfunclist{strcmpi(ntype,ntypelist)}(imsize,corr,cfreq);

% get required noise image by IFFT
I = real(ifft2(ifftshift(F)));

% map image elements to [0,1] then convert to uint8
I = (I - min(I(:))) / range(I(:));

% binary process
if any(validatestring(ntype,bproclist))
    I = I > median(I(:));
else
    I = im2uint8(I);
end

end

function mask = white(varargin)
% white mask mean every pixel is independent, there is no dependence at
% all. So, the spectrum should be equal at every position.
mask = 1;
end

function mask = pinklike(sz,corr,varargin)
% pinklike noise is a set of noise that have spectrum as power of 1/f. The
% original pink noise is 1/f, with the power term getting higher the
% dependency between pixels are getting stronger.

% generate coordinates matrix and frequency map
switch numel(sz)
    case 1
        [c,r] = meshgrid(1:sz,1:sz);
        freq = sqrt((c-(sz+1)/2).^2 + (r-(sz+1)/2).^2);
    case 2
        [c,r] = meshgrid(1:sz(1),1:sz(2));
        freq  = sqrt((c-(sz(1)+1)/2).^2 + (r-(sz(2)+1)/2).^2);
    otherwise
        error('Input argument SZ should have at most 2 elements.');
end

% generate mask according to power law (add EPS to avoid Inf)
mask = power(freq+eps,-corr);

end

function mask = cow(sz,corr,cfreq, varargin)
% cow noise is a band-pass noise. when binarize the noise image, a cow
% texture would show.

% generate coordinates matrix and frequency map
switch numel(sz)
    case 1
        [c,r] = meshgrid(1:sz,1:sz);
        freq = sqrt((c-(sz+1)/2).^2 + (r-(sz+1)/2).^2);
    case 2
        [c,r] = meshgrid(1:sz(1),1:sz(2));
        freq  = sqrt((c-(sz(1)+1)/2).^2 + (r-(sz(2)+1)/2).^2);
    otherwise
        error('Input argument SZ should have at most 2 elements.');
end

% generate mask according to power law (add EPS to avoid Inf)
mask = exp(-power((freq-cfreq)*corr,2));

end
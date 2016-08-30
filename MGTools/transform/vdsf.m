function v = vdsf(v,fsz,nfrm)
% VDSF downsample video by filtering the video with spacial-frequncy domain
% filter on spacial axes and 1D Gaussian filter on temporal axes.
%
%   V = VDSF(V,FSZ,NFRM) downsample video V to a video with frame size FSZ
%   and frame number NFRM.
%
% see also vds, movresize, sample, and downsample.
%
% MooGu Z. <hzhu@case.edu>
% Jun 5, 2014 - Version 0.1

% Rate between filter size and its standard deviation
filterSizeFactor = 3; 

% Size of Original Video
[orsz,ocsz,onfrm] = size(v);
% Frame Size of Downsample Video
if numel(fsz)==1, rsz = fsz; csz = fsz;
else rsz = fsz(1); csz = fsz(2); end
% Fix the shift problem
if (orsz - rsz > 2) && (ocsz - csz > 2)
    swShift = true;
    rsz = rsz + 2;
    csz = csz + 2;
else
    swShift = false;
end

assert(rem(orsz-rsz,2)==0 && rem(ocsz-csz,2)==0,...
    'Resolution of downsample should have same parity of original one.');

% Cutoff Power of Gaussian Filter
tCutoff = 1e-5;
% Estimate Standard Deviation of Gaussian Filter
tsigma = norminv(1-tCutoff/2,0,1) / (nfrm*pi); % Time

% Calculate compensate ration in FFT transform
cRatio = (orsz*ocsz) / (rsz*csz);
% Calculate center area coordinates
rArea = (orsz-rsz)/2 + 1 : orsz - (orsz-rsz)/2;
cArea = (ocsz-csz)/2 + 1 : ocsz - (ocsz-csz)/2;
% Generate Mask
mask = sFilterMask([orsz,ocsz],.2*(rsz+csz));
% filtering on spacial axes
buffer = zeros(rsz,csz,onfrm);
for t = 1 : onfrm
    masked = mask .* fftshift(fft2(v(:,:,t)));
    buffer(:,:,t) = real(ifft2(ifftshift(masked(rArea,cArea)))) / cRatio;
    % buffer(:,:,t) = real(ifft2(ifftshift(masked)));
end

v = buffer;

% filtering on temporal axes
if onfrm == nfrm 
    % do nothing
elseif mod(onfrm,nfrm) == 0
    filterSize = floor(filterSizeFactor * tsigma * onfrm);
    f = tFilterFunc((-filterSize:filterSize)/onfrm,tsigma);
    normFactor = conv(ones(onfrm,1),f,'same');
    for r = 1 : rsz
        for c = 1 : csz
            v(r,c,:) = conv(reshape(v(r,c,:),onfrm,1),f,'same') ...
                ./ normFactor;
        end
    end
    sampleIndex = 1 : onfrm/nfrm : onfrm;
    v = v(:,:,sampleIndex);
else
    buffer = zeros(rsz,csz,nfrm);
    opos = linspace(0,1,onfrm);  % Original Pixel Position
    dpos = linspace(0,1,nfrm);   % Downsample Pixel Position
    for t = 1 : nfrm
        f = reshape(tFilterFunc(opos-dpos(t),tsigma),1,1,onfrm);
        for r = 1 : rsz
            for c = 1 : csz
                buffer(r,c,t) = sum(v(r,c,:).*f);
            end
        end
    end
    v = buffer;
end

if swShift
    v = v(2:end-1,2:end-1,:);
end

end

function mask = sFilterMask(sz,fstop)
% SFILTERMASK create a mask with size SZ and stop frequency FSTOP. Here we
% use 4-th order low-pass filter utilized in Olshausen's Nature paper

assert(length(sz)==2,'Input argument SZ should have 2 elements.');

% Generate coordinates of each point
[c,r] = meshgrid(1:sz(1),1:sz(2));
% Calculate mask factor on each coordinate
freq = sqrt((r-(sz(1)+1)/2).^2 + (c-(sz(2)+1)/2).^2);
mask = exp(-(freq/fstop).^4);

end

function f = tFilterFunc(x,sigma)
% FILTERFUNCTION map numbers to filter function value

f = exp(-(x.^2)/(2*sigma^2));   % Gaussian Distribution PDF
f = f / sum(f(:));              % Normalization

end
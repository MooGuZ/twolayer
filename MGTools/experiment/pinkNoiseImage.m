function img = pinkNoiseImage(sz)
% PINKNOISEIMAGE generate a image with expected size SZ containing noise
% with 'pink' power spectrum.
%
% MooGu Z. <hzhu@case.edu>
% Jan 15, 2015 - Version 0.1

% Generate white noise image with the expected size
rng('shuffle');
img = randn(sz);

% Calculate the spacial frequency spectrum of the white noise image
spectrum = fftshift(fft2(img));

% Apply 1/f power spectrum rule with masking
spectrum = sFilterMask(sz) .* spectrum;

% Get pink noise image by IFFT
img = real(ifft2(ifftshift(spectrum)));

% map image elements to [0,1]
img = im2uint8(img - min(img(:))) / (max(img(:)) - min(img(:)));

end

function mask = sFilterMask(sz)
% SFILTERMASK create a mask with size SZ and stop frequency FSTOP. Here we
% use 4-th order low-pass filter utilized in Olshausen's Nature paper

% Check the input argument and modify it, if necessary
if length(sz) == 1
    sz = [sz,sz];
elseif length(sz) > 2
    error('Input argument SZ should have at most 2 elements.');
end
% Generate coordinates matrix
[c,r] = meshgrid(1:sz(1),1:sz(2));
% Calculate frequency of each postion
freq = sqrt((r-(sz(1)+1)/2).^2 + (c-(sz(2)+1)/2).^2);
% Generate mask according to power law (add EPS to avoid devided by zero)
mask = 1 ./ sqrt(freq+eps);

end
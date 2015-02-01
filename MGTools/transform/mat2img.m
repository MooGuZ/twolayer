function I = mat2img(M,anchor)
% MAT2IMG transform matrix into colorful image. This transformation  
% modulate (saturation,brightness) by the norm of numbers in the matrix
% and hue modulated by phase of numbers in the matrix. In this way, complex
% number would map to a cone surface in HSV color space, while real number
% matrix corresponds to a [cyan,black,red] color region.
%
% Usage:
% 
% I = MAT2IMG(M,ANCHOR), where M is the 2D matrix and ANCHOR means the
% value of (S,V) pair that would be mapped from 1.0 (if M is complex, 
% ANCHOR corresponds to the complex number with norm of 1).
%
% MooGu Z. <hzhu@case.edu>
% 2015.01.23 - Version 0.1

assert(numel(size(M))==2,'Matrix has and only has two dimensions.');
% Set default CPOINT
if ~exist('anchor','var'), anchor = .9; end
% Check input argument
assert((anchor > 0) && (anchor < 1), ...
    'The anchor should be in the range (0,1).');
% Generate Image with same size as M
I = zeros([size(M),3]);
% Calculate Hue value of each element
if isreal(M)
    I(M(:)<0) = .5;
else
    I(:,:,1)  = (wrapToPi(angle(M)-pi)/pi+1) / 2;
end
% Calculate (Saturation,Value) pairs
interMap = atan(tan(anchor*(pi/2))*abs(M))*(2/pi);
% Tweak the curve
area = (abs(M) <= 1);
interMap(area) = anchor * (abs(M(area))).^(2/3);
% Fill up (S,V) pairs
I(:,:,2:3) = repmat(interMap,[1,1,2]);
% Convert Image into RGB Space
I = hsv2rgb(I);

end


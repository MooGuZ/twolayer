function I = mat2img(M,method,th)
% MAT2IMG transform matrix into colorful image. This transformation  
% modulate (saturation,brightness) by the norm of numbers in the matrix
% and hue modulated by phase of numbers in the matrix. In this way, complex
% number would map to a cone surface in HSV color space, while real number
% matrix corresponds to a [cyan,black,red] color region.
%
% Usage:
% 
% I = MAT2IMG(M,[METHOD,TH]), where M is the 2D matrix and METHOD is
% mapping method utilized in this transformation process. There are two
% method available now : ATAN and SIGMOID. If SIGMOID is active, threshold
% TH can be specified.
%
% MooGu Z. <hzhu@case.edu>
% 2015.01.23 - Version 0.1

assert(numel(size(M))==2,'Matrix shoule has and only has two dimensions.');
% Set default CPOINT
if ~exist('method','var'), method = 'atan'; end
if ~exist('th','var'), th = .3; end
% Generate Image with same size as M
I = zeros([size(M),3]);
% Setup (Saturation,Value) pairs
switch lower(method)
    case {'atan'}
        % Calculate (Saturation,Value) pairs
        interMap = atan(tan(.9*(pi/2))*abs(M))*(2/pi);
        % Tweak the curve
        area = (abs(M) <= 1);
        interMap(area) = .9 * (abs(M(area))).^(2/3);
        % Fill up (S,V) pairs
        I(:,:,2:3) = repmat(interMap,[1,1,2]);
    case {'sigmoid'}
        I(:,:,2:3) = repmat(1./(1+exp(-(abs(M)-th)*10)),[1,1,2]);
    otherwise
        error('the method is undefined!');
end
% Setup Hue Value
if isreal(M)
    I(M(:)<0) = .5;
else
    I(:,:,1)  = (wrapToPi(angle(M)-pi)/pi+1) / 2;
end
% Convert Image into RGB Space
I = hsv2rgb(I);

end


function anim2gif(frames,sz,fname,swNorm,frate)
% ANIM2GIF would convert animation stored in column frames into gif file
%
%   anim2gif(frames,sz,fname,swNorm)
%
% MooGu Z. <hzhu@case.edu>
% April 22, 2014 - Version 0.2

if exist('frate','var')
    delay = 1/frate;
else
    delay = 0.0417;
end

[npixels,nframe] = size(frames);

if ~exist('sz','var')
    sz = [floor(sqrt(npixels)),floor(sqrt(npixels))];
    assert(prod(sz) == npixels,...
        '[ANIM2GIF] This function need square images or specified frame size');
end

if ~exist('fname','var')
    fname = ['./fig/',datestr(now),'.gif'];
end

% Normalize
if exist('swNorm','var') && swNorm
    Min = min(frames(:));
    Max = max(frames(:));
    frames = (frames - Min) / (Max - Min);
end

[I,map] = gray2ind(reshape(frames(:,1),sz(1),sz(2)),256);
imwrite(I,map,fname,'gif','DelayTime',delay);

for i = 2 : nframe
    [I,map] = gray2ind(reshape(frames(:,i),sz(1),sz(2)),256);
    imwrite(I,map,fname,'gif','WriteMode','append','DelayTime',delay);
end

end
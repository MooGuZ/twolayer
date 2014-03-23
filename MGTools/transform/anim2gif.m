function anim2gif(frames,sz)
% ANIM2GIF would convert animation stored in column frames into gif file
% under directory of 'fig'

[npixels,nframe] = size(frames);

if ~exist('sz','var')
    sz = [floor(sqrt(npixels)),floor(sqrt(npixels))];
    assert(prod(sz) == npixels,...
        '[FRM2GIF] This function need square images or specified frame size');
end

fname = ['./fig/',datestr(now),'.gif'];

% Normalize
frames = bsxfun(@minus,frames,min(frames));
frames = bsxfun(@rdivide,frames,max(frames));

[I,map] = gray2ind(reshape(frames(:,1),sz(1),sz(2)),256);
imwrite(I,map,fname,'gif','DelayTime',0.0417);

for i = 2 : nframe
    [I,map] = gray2ind(reshape(frames(:,i),sz(1),sz(2)),256);
    imwrite(I,map,fname,'gif','WriteMode','append','DelayTime',0.0417);
end

end
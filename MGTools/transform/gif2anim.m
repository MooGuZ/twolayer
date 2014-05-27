function frames = gif2anim(fname)
% GIF2ANIM read GIF file into animation in column frames

[I,cmap] = imread(fname,'Frames','all');

isz = size(I);
assert(isz(3)==1,'[GIF2ANIM] GIF file cannot be recognized!');

% Convert Color Map to YCbCr Space
cmap = rgb2ycbcr(cmap);
cmap = cmap(:,1);

% Reshape Index Matrix
I = reshape(I,isz(1)*isz(2),isz(4));
% Construct Frames
frames = cmap(double(I)+1);

end

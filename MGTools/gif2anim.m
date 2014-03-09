function frames = gif2anim(fname)
% GIF2ANIM read GIF file into animation in column frames

[I,map] = imread(fname,'Frames','all');

isz = size(I);
frames = zeros(isz(1)*isz(2),isz(4));

for i = 1 : isz(4)
    img = double(I(:,:,1,i)) + 1;
    ycbcr = rgb2ycbcr(map(img(:),:));
    frames(:,i) = ycbcr(:,1);
end

end

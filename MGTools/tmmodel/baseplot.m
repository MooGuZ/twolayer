function I = baseplot(pbases,tbases,bsz)
% BASEPLOT draw image of base

% dimensions of bases
[npixel,npat] = size(pbases);
[~,ntrans]    = size(tbases);
% set default value
if ~exist('bsz','var'), bsz = round(sqrt(npixel)) * [1,1]; end
% check input arguments
assert(size(tbases,1) == npixel, 'Bases are not match!');
assert(prod(bsz) == npixel, 'Bases size need to be specified!');

% interval size between bases
bisz = 2;
% interval size between sections
sisz = 4;
% calculate image size and initialize it
M = zeros((npat + 1) * (bsz(1) + bisz) + 3*sisz - 2*bisz, ...
    (ntrans + 1) * (bsz(2) + bisz) + 3*sisz - 2*bisz);
% define coordinates function
up    = @(i) (bsz(1) + bisz) * i + 2*sisz - bisz + 1;
down  = @(i) (bsz(1) + bisz) * (i + 1) + 2*(sisz - bisz);
left  = @(j) (bsz(2) + bisz) * j + 2*sisz - bisz + 1;
right = @(j) (bsz(2) + bisz) * (j + 1) + 2*(sisz - bisz);

% draw pattern bases
for i = 1 : npat
    M(up(i):down(i),sisz+1:sisz+bsz(2)) = reshape(pbases(:,i),bsz);
end
% draw transform bases
for j = 1 : ntrans
    M(sisz+1:sisz+bsz(1),left(j):right(j)) = ...
        exp(1j*reshape(tbases(:,j),bsz));
end
% draw complex bases
for i = 1 : npat
    for j = 1 : ntrans
        M(up(i):down(i),left(j):right(j)) = ...
            reshape(pbases(:,i),bsz) .* exp(1j*reshape(tbases(:,j),bsz));
    end
end

% convert matrix to image
I  = mat2img(M);

end
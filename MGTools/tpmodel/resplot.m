function I = resplot(R)
% RESPLOT draw a image show the respond of TP-Model

% check input arguments
assert(isnumeric(R) && ~isreal(R) && numel(size(R)) == 3, ...
    'Input is 3-D complex matrix!');

% block size for each respond
bsz = 2;
% interval between responds
bisz = 1;
% interval between frames
fisz = 4;

% check size of responds
[npat,ntrans,nfrm] = size(R);
% define arrangement of frames
ncol = ceil(sqrt(nfrm*npat/ntrans));
nrow = ceil(nfrm/ncol);
ncol = ceil(nfrm/nrow);
% calculate size of plot and initilize it
I = ones(nrow * ((npat * (bsz + bisz) + bisz) + fisz) + fisz, ...
    ncol * ((ntrans * (bsz + bisz) + bisz) + fisz) + fisz, 3);
% define coordinate function for each frame
left  = @(f) rem(f-1,ncol) * ((ntrans * (bsz + bisz) + bisz) + fisz) + fisz + 1;
right = @(f) (rem(f-1,ncol) + 1) * ((ntrans * (bsz + bisz) + bisz) + fisz);
up    = @(f) (ceil(f/ncol) - 1) * ((npat * (bsz + bisz) + bisz) + fisz) + fisz + 1;
down  = @(f) ceil(f/ncol) * ((npat * (bsz + bisz) + bisz) + fisz);
% initilize plot buffer for each frame
buffer = zeros(npat * (bsz + bisz) + bisz,ntrans * (bsz + bisz) + bisz);
% generate kron kernel
kern = [ones(bsz,1);zeros(bisz,1)];
kern = kern * kern';
% draw each frame
for f = 1 : nfrm
    % fill responds to buffer
    buffer(bisz+1:end,bisz+1:end) = kron(R(:,:,f),kern);
    % plot current respond frame
%     I(up(f):down(f),left(f):right(f),:) = mat2img(buffer);
    I(up(f):down(f),left(f):right(f),:) = mat2img(abs(buffer));
end

end
function I = resplot(beta,theta,bia)
% RESPLOT draw a image show the respond of TP-Model

% block size for each respond
bsz = 2;
% interval between responds
bisz = 1;
% interval between frames
fisz = 4;

% check size of responds
[npat,ntrans,nfrm] = size(beta);
% check dimensional consistancy
assert(all([npat,ntrans,nfrm] == size(theta)), ...
    'dimension of THETA does not match!');
assert(all([npat,nfrm] == size(bia)),'dimension of BIA does not match!');
% define arrangement of frames
ncol = ceil(sqrt(nfrm*npat/(ntrans+1)));
nrow = ceil(nfrm/ncol);
ncol = ceil(nfrm/nrow);
% calculate size of plot and initilize it
I = ones(nrow * ((npat * (bsz + bisz) + bisz) + fisz) + fisz, ...
    ncol * (((ntrans + 1) * (bsz + bisz) + bisz) + fisz) + fisz, 3);
% define coordinate function for each frame
left  = @(f) rem(f-1,ncol) * ...
    (((ntrans + 1) * (bsz + bisz) + bisz) + fisz) + fisz + 1;
right = @(f) (rem(f-1,ncol) + 1) * ...
    (((ntrans + 1) * (bsz + bisz) + bisz) + fisz);
up    = @(f) (ceil(f/ncol) - 1) * ...
    ((npat * (bsz + bisz) + bisz) + fisz) + fisz + 1;
down  = @(f) ceil(f/ncol) * ...
    ((npat * (bsz + bisz) + bisz) + fisz);
% initilize plot buffer for each frame
buffer = zeros(npat * (bsz + bisz) + bisz,(ntrans + 1) * ...
    (bsz + bisz) + bisz);
% generate kron kernel
kern = [ones(bsz,1);zeros(bisz,1)];
kern = kern * kern';
% assemble responds
R = [reshape(bia,[npat,1,nfrm]),beta.*exp(1j*theta)];
% draw each frame
for f = 1 : nfrm
    % fill responds to buffer
    buffer(bisz+1:end,bisz+1:end) = kron(R(:,:,f),kern);
    % plot current respond frame
    I(up(f):down(f),left(f):right(f),:) = mat2img(buffer);
end

% if this function is running alone, just show the image and return the
% figure handle.
if nargout == 0
    fig = figure();
    imshow(I);
    % return the figure handle instead of image matrix
    I = fig;
end

end

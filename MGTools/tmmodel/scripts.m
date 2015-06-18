% This file contains all temporal scripts for tmmodel

ename = 'motion-hmsparse';

%% save reconstruction
anim2gif(video.v,sprintf('~/Desktop/fig/%s-org.gif',ename));
anim2gif(rec,sprintf('~/Desktop/fig/%s-rec.gif',ename));

%% draw mask bases
nrow = 6;
ncol = 11;
res  = [32,32];

bisz = 3; % number of pixel between bases

I = ones(bisz + nrow * (bisz + res(1)), bisz + ncol * (bisz + res(2)),3);

% define coordinate function for each frame
left  = @(f) rem(f-1,ncol) * (res(1) + bisz) + bisz + 1;
right = @(f) (rem(f-1,ncol) + 1) * (res(1) + bisz);
up    = @(f) (ceil(f/ncol) - 1) * (res(2) + bisz) + bisz + 1;
down  = @(f) ceil(f/ncol) * (res(2) + bisz);

for i = 1 : size(m.alpha, 2)
    I(up(i):down(i),left(i):right(i),:) = mat2img(reshape(m.alpha(:,i),res));
end

imwrite(I,sprintf('~/Desktop/fig/%s-mask.png',ename),'png');

%% draw transformation bases
nrow = 4;
ncol = 4;
res  = [32,32];

bisz = 3; % number of pixel between bases

I = ones(bisz + nrow * (bisz + res(1)), bisz + ncol * (bisz + res(2)),3);

% define coordinate function for each frame
left  = @(f) rem(f-1,ncol) * (res(1) + bisz) + bisz + 1;
right = @(f) (rem(f-1,ncol) + 1) * (res(1) + bisz);
up    = @(f) (ceil(f/ncol) - 1) * (res(2) + bisz) + bisz + 1;
down  = @(f) ceil(f/ncol) * (res(2) + bisz);

for i = 1 : size(m.phi, 2)
    I(up(i):down(i),left(i):right(i),:) = ...
        mat2img(exp(1j*reshape(m.phi(:,i),res)));
end

imwrite(I,sprintf('~/Desktop/fig/%s-trans.png',ename),'png');

%% save beta distribution
hist(m.beta(:),1000);
print(gcf,'-dpng',sprintf('~/Desktop/fig/%s-beta-dist.png',ename));

%% save theta distribution
hist(m.theta(:),1000);
print(gcf,'-dpng',sprintf('~/Desktop/fig/%s-theta-dist.png',ename));

%% save sample animation and responds
ind = 25 : 48;
sname = 'sampleA';
anim2gif(rec(:,ind),sprintf('~/Desktop/fig/%s-%s.gif',ename,sname));
imwrite(resplot(m.beta(:,:,ind)*10, m.theta(:,:,ind), m.bia(:,ind)), ...
    sprintf('~/Desktop/fig/%s-%s-resp.png',ename,sname),'png');

% END
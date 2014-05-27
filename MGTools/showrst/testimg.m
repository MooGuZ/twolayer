function img = testimg(sz)
% Generate Test Image for Simulation Motion Code learned by C&O model

alpha = 6;

[x,y] = meshgrid(1:sz,1:sz);
% Centralization
x = bsxfun(@minus,x,mean(x,2));
y = bsxfun(@minus,y,mean(y));

% Generating Test Image
R2   = x.^2 + y.^2;
Rmax = max(R2(:));
img  = .5 + .5*sin(cos(alpha*pi*R2/Rmax)*pi);

end
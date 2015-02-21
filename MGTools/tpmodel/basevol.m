function [anim,fres] = basevol(video,model)
% BASEEVOL create animation shows the evolution of bases in current status

% set number of animation frames
nfrm = 24;
% set evolution step size
step = 1;

% get dimensionality of data
[npixel,nframe] = size(video.v);
npattern = size(model.alpha,2);
ntrans   = size(model.phi,2);
% convert model to parameters
alpha = reshape(model.alpha,[npixel,npattern,1,1]);
phi   = reshape(model.phi,[npixel,1,ntrans,1]);
beta  = reshape(model.beta,[1,npattern,ntrans,nframe]);
theta = reshape(model.theta,[1,npattern,ntrans,nframe]);
bia   = reshape(model.bia,[1,npattern,1,nframe]);
v     = reshape(video.v,[npixel,1,1,nframe]);
% get sigma setting from model
sigma = model.sigma;

% calculate error in current status
delta = v - genmodel(alpha,phi,beta,theta,bia);
% calculate derivatives of bases
[dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia,delta, ...
        sigma,video.res,true,true);
    
% generate first frame
frame = baseplot(reshape(alpha,[npixel,npattern]),reshape(phi,[npixel,ntrans]));
% get dimension of first frame
[height,width,~] = size(frame);
% initilize animation
anim = zeros(height*width,nfrm,3);
anim(:,1,:) = reshape(frame,[height*width,1,3]);
% generate following frames
for i = 2 : nfrm
    alpha = alpha - step * dAlpha;
    phi   = phi - step * dPhi;
    anim(:,i,:) = reshape(baseplot(reshape(alpha,[npixel,npattern]), ...
        reshape(phi,[npixel,ntrans])),[height*width,1,3]);
end

% generate returen value
fres = [height,width];

end

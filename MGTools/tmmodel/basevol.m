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
[alpha,phi,beta,theta,bia] = m2p(model);
% reshape video for calculating
v = reshape(video.v,[npixel,1,1,nframe]);
% get sigma setting from model
sigma = model.sigma;
% get control parameters from model
ctrl  = model.ctrl;

% calculate error in current status
delta = v - genmodel(alpha,phi,beta,theta,bia);
% calculate derivatives of bases
[dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia,delta, ...
        sigma,ctrl,video.res);
    
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

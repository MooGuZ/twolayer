function [alpha,phi,beta,theta,bia] = m2p(model)
% M2P convert model to parameter list by reshaping bases and responds to
% corresponding 4-D matrix

% get dimension info
[npixel,npattern] = size(model.alpha);
ntrans = size(model.phi,2);
% reshape bases for calculation convenience
alpha = reshape(model.alpha,[npixel,npattern,1,1]);
phi   = reshape(wrapToPi(model.phi),[npixel,1,ntrans,1]);
% check existence of responds
if isfield(model,'beta')
    nframe = size(model.beta,3);
    beta   = reshape(model.beta,[1,npattern,ntrans,nframe]);
    theta  = reshape(wrapToPi(model.theta),[1,npattern,ntrans,nframe]);
    bia    = reshape(model.bia,[1,npattern,1,nframe]);
else
    beta   = [];
    theta  = [];
    bia    = [];
end

end
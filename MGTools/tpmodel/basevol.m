function [anim,res] = basevol(v,model)
% BASEEVOL create animation shows the evolution of bases in current status

% set number of animation frames
nfrm = 24;
% set evolution step size
step = 1;

% convert model to parameters
[alpha,phi,beta,theta,bia] = tpmodel2data(model);
% get sigma setting from model
sigma = model.sigma;
% get dimensions
[npixel,~,npat] = size(alpha);
[~,~,~,ntrans]  = size(beta);
% estimate resolution
res = round(sqrt(npixel)) * [1,1];

% calculate error in current status
delta = v - genModel(alpha,phi,beta,theta,bia);
% calculate derivatives of bases
[dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia,delta, ...
        sigma,res,true,true);
    
% generate first frame
frame = baseplot(reshape(alpha,[npixel,npat]),reshape(phi,[npixel,ntrans]));
% get dimension of first frame
[height,width,~] = size(frame);
% initilize animation
anim = zeros(height*width,nfrm,3);
anim(:,1,:) = reshape(frame,[height*width,1,3]);
% generate following frames
for i = 2 : nfrm
    alpha = alpha + step * dAlpha;
    phi   = phi + step * dPhi;
    anim(:,i,:) = reshape(baseplot(reshape(alpha,[npixel,npat]), ...
        reshape(phi,[npixel,ntrans])),[height*width,1,3]);
end

% generate returen value
res = [height,width];

end

function [alpha,phi,beta,theta,bia] = tpmodel2data(Model)
% Scale of each dimension
[npixel,ntrans] = size(Model.transBase);
[~,npattern]    = size(Model.patBase);
% Transform Data into 4D
phi    = reshape(Model.transBase,[npixel,1,1,ntrans]);
alpha  = reshape(Model.patBase,[npixel,1,npattern,1]);
theta  = permute(Model.transCoefficient,[4,3,1,2]);
beta   = permute(Model.patCoefficient,[4,3,1,2]);
bia    = permute(Model.bia,[4,2,1,3]);
end

function obj = objFunc(alpha,phi,beta,theta,~,delta,sigma,ffindex,res)
npixel   = size(alpha,1);
nframe   = size(beta,2);
npattern = size(theta,3);
ntrans   = size(phi,4);
nrow     = res(1);
ncol     = res(2);
% Difference of each segment along time axis of theta
segDiff = wrapToPi(diff(theta,1,2));
segDiff(:,ffindex(2:end)-1,:,:) = 0;
% Reshape Bases for the convenience of calculation
alpha = reshape(alpha,[nrow,ncol,npattern]);
phi   = reshape(phi,[nrow,ncol,ntrans]);
% Object Values
obj.noise   = sum(delta(:).^2) / (sigma.noise^2 * npixel * nframe);
obj.sparse  = sum(log(1+(beta(:)/sigma.sparse).^2)) / ...
    (nframe * npattern * ntrans);
obj.slow    = sum(segDiff(:).^2) / ...
    (sigma.slow^2 * nframe * npattern * ntrans);
obj.smpat   = ...
    (sum(reshape(diff(alpha,1,1).^2,[(npixel-ncol)*npattern,1])) ...
    + sum(reshape(diff(alpha,1,2).^2,[(npixel-nrow)*npattern,1]))) ...
    / (2 * sigma.smpat^2 * npixel * npattern);
obj.smtrans = ...
    (sum(reshape(wrapToPi(diff(phi,1,1)).^2,[(npixel-ncol)*ntrans,1])) ...
    + sum(reshape(wrapToPi(diff(phi,1,2)).^2,[(npixel-nrow)*ntrans,1]))) ...
    / (2 * sigma.smtrans^2 * npixel * ntrans);
obj.value   = obj.noise + obj.sparse + obj.slow + ...
    obj.smpat + obj.smtrans;
end

function v = genModel(alpha,phi,beta,theta,bia)
v = sum(bsxfun(@times,alpha,bsxfun(@plus,bia, ...
        sum(bsxfun(@times,beta,cos(bsxfun(@minus,phi,theta))),4))),3);
end

function [dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia,delta, ...
    sigma,res,swAlpha,swPhi)
if ~(swAlpha || swPhi), return; end
% Initialize derivatives
dAlpha = 0;
dPhi   = 0;
% Common term in derivatives of alpha and phi
phase  = bsxfun(@minus,phi,theta);
% Noise part of derivatives of alpha and phi
if swAlpha
    dAlpha = -sum(bsxfun(@times,delta,bsxfun(@plus,bia, ...
        sum(bsxfun(@times,beta,cos(phase)),4))),2) / sigma.noise^2;
end
if swPhi
    dPhi   = sum(bsxfun(@times,delta, ...
        sum(bsxfun(@times,alpha,beta).*sin(phase),3)),2) / sigma.noise^2;
end
% Reshape alpha and phi to 3D matrix for calculation convenience
if swAlpha, alpha = reshape(alpha,[res,size(beta,3)]); end
if swPhi,   phi   = reshape(phi,[res,size(theta,4)]);  end
% Calculate smoothness part derivatives
if swAlpha
    dAlpha = dAlpha - reshape(diff(padarray(alpha,[1,0,0],'replicate'),2,1) ...
        + diff(padarray(alpha,[0,1,0],'replicate'),2,2),size(dAlpha)) ...
        / (2 * numel(alpha) * sigma.smpat^2 / numel(delta));
end
if swPhi
    dPhi = dPhi - reshape(diff(padarray(phi,[1,0,0],'replicate'),2,1) ...
        + diff(padarray(phi,[0,1,0],'replicate'),2,2),size(dPhi)) ...
        / (2 * numel(phi) * sigma.smtrans^2 / numel(delta));
end
% Normalize Gradients
normFactor = sqrt(sum(dAlpha(:).^2) + sum(dAlpha(:).^2));
if swAlpha, dAlpha = dAlpha / normFactor; end
if swPhi,   dPhi   = dPhi   / normFactor; end
end
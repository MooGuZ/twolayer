function [dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia,delta, ...
    sigma,res,swAlpha,swPhi)
if ~(swAlpha || swPhi), return; end
% Initialize derivatives
if gpuDeviceCount == 0
    dAlpha = zeros(size(alpha));
    dPhi   = zeros(size(phi));
else
    dAlpha = zeros(size(alpha),'gsingle');
    dPhi   = zeros(size(phi),'gsingle');
end
% calculate noise part of derivatives frame by frame
for f = 1 : size(beta,4)
    % Common term in derivatives of alpha and phi
    phase  = bsxfun(@minus,phi,theta(:,:,:,f));
    % Noise part of derivatives of alpha and phi
    if swAlpha
        dAlpha = dAlpha - bsxfun(@times,delta(:,:,:,f),bsxfun(@plus, ...
            bia(:,:,:,f),sum(bsxfun(@times,beta(:,:,:,f),cos(phase)),3))) ;
    end
    if swPhi
        dPhi = dPhi + bsxfun(@times,delta(:,:,:,f), ...
            sum(bsxfun(@times,alpha,beta(:,:,:,f)).*sin(phase),2));
    end
end
% Reshape alpha and phi to 3D matrix for calculation convenience
if swAlpha, alpha = reshape(alpha,[res,size(beta,2)]); end
if swPhi,   phi   = reshape(phi,[res,size(theta,3)]);  end
% Calculate smoothness part derivatives
if swAlpha
    dAlpha = dAlpha / sigma.noise^2 - reshape(diff(padarray(alpha,[1,0,0],'replicate'),2,1) ...
        + diff(padarray(alpha,[0,1,0],'replicate'),2,2),size(dAlpha)) ...
        / (2 * numel(alpha) * sigma.smpat^2 / numel(delta));
end
if swPhi
    dPhi = dPhi / sigma.noise^2 - reshape( ...
        diff(wrapToPi(diff(padarray(phi,[1,0,0],'replicate'),1,1)),1,1) ...
        + diff(wrapToPi(diff(padarray(phi,[0,1,0],'replicate'),1,2)),1,2),size(dPhi)) ...
        / (2 * numel(phi) * sigma.smtrans^2 / numel(delta));
end
% Normalize Gradients
normFactor = sqrt(sum(dAlpha(:).^2) + sum(dPhi(:).^2));
if swAlpha, dAlpha = dAlpha / normFactor; end
if swPhi,   dPhi   = dPhi   / normFactor; end
end

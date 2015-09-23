function [dAlpha,dPhi,normFactor] = ...
    dBase(alpha,phi,beta,theta,bia,delta,sigma,ctrl,res)

% short path for easy case
if ~(ctrl.swPatOpt || ctrl.swTransOpt)
	dAlpha = 0;
	dPhi = 0;
	normFactor = 1;
	return;
end

% derivative of alpha
if ctrl.swPatOpt
    % initialize derivatives of noise part
    dNoise = zeros(size(alpha));
    if gpuDeviceCount ~= 0
        dNoise = gsingle(dNoise);
    end
    % calculate noise part of derivatives of frame by frame
    for f = 1 : size(beta,4)
        phase  = bsxfun(@minus,phi,theta(:,:,:,f));
        % original derivatives
        dNoise = dNoise - bsxfun(@times,delta(:,:,:,f),bsxfun(@plus, ...
            bia(:,:,:,f),sum(bsxfun(@times,beta(:,:,:,f),cos(phase)),3)));
    end
    % modify derivative for alpha-normalization
    if ctrl.swANorm && ctrl.swNegCut
        dNoise = compose(repmat(-sum(dNoise .* alpha,1) / ctrl.anorm, ...
            [prod(res),1]),dNoise,alpha);
    end
    % reshape alpha to calculate smoothness part derivatives
    alpha = reshape(alpha,[res,size(beta,2)]);
    % intermediate results
    drow = diff(alpha,1,1);
    dcol = diff(alpha,1,2);
    % calculate original smooth derivatives
    dSmooth = - (diff(padarray(drow,[1,0]),1,1) ...
        + diff(padarray(dcol,[0,1]),1,2));
    % modify derivative for alpha-normalization
    if ctrl.swANorm && ctrl.swNegCut
        dSmooth = compose(repmat(-(sum(sum(drow.^2,1),2) ...
            + sum(sum(dcol.^2,1),2)) / ctrl.anorm,[res,1]),dSmooth,alpha);
    end
    % reshape alpha back
    alpha = reshape(alpha,[prod(res),size(beta,2)]);
    % assemble final derivatives of alpha
    dAlpha = dNoise / sigma.noise^2 + reshape(dSmooth,size(alpha)) ...
        / (2 * numel(alpha) * sigma.smpat^2 / numel(delta));
else
    dAlpha = 0;
end

% derivative of phi
if ctrl.swTransOpt
    % initialize derivatives of noise part
    dNoise = zeros(size(phi));
    if gpuDeviceCount ~= 0
        dNoise = gsingle(dNoise);
    end
    % calculate noise part of derivatives of frame by frame
    for f = 1 : size(beta,4)
        phase  = bsxfun(@minus,phi,theta(:,:,:,f));
        dNoise = dNoise + bsxfun(@times,delta(:,:,:,f), ...
            sum(bsxfun(@times,alpha,beta(:,:,:,f)).*sin(phase),2));
    end
    % reshape phi to calculate derivatives of smoothness
    phi  = reshape(phi,[res,size(theta,3)]);
    % intermediate results
    drow = wrapToPi(diff(phi,1,1));
    dcol = wrapToPi(diff(phi,1,2));
    % calculate derivative of smoothness
    dSmooth = - (diff(padarray(drow,[1,0]),1,1) ...
        + diff(padarray(dcol,[0,1]),1,2));
    % reshape alpha back
    phi = reshape(phi,[prod(res),1,size(theta,3)]);
    % assemble final derivatives of phi
    dPhi = dNoise / sigma.noise^2 + reshape(dSmooth,size(phi)) ...
        / (2 * numel(phi) * sigma.smtrans^2 / numel(delta));
else
    dPhi = 0;
end

% % Initialize derivatives
% if gpuDeviceCount == 0
%     dAlpha   = zeros(size(alpha));
%     dPhi     = zeros(size(phi));
% else
%     dAlpha   = zeros(size(alpha),'gsingle');
%     dPhi     = zeros(size(phi),'gsingle');
% end
% 
% % calculate noise part of derivatives of frame by frame (alpha normalized)
% for f = 1 : size(beta,4)
%     % Common term in derivatives of alpha and phi
%     phase  = bsxfun(@minus,phi,theta(:,:,:,f));
%     % Noise part of derivatives of alpha and phi
%     if ctrl.swPatOpt
%         % original derivatives without normalization
%         dAlpha = dAlpha - bsxfun(@times,delta(:,:,:,f),bsxfun(@plus, ...
%             bia(:,:,:,f),sum(bsxfun(@times,beta(:,:,:,f),cos(phase)),3)));
%     end
%     if ctrl.swTransOpt
%         dPhi = dPhi + bsxfun(@times,delta(:,:,:,f), ...
%             sum(bsxfun(@times,alpha,beta(:,:,:,f)).*sin(phase),2));
%     end
% end
% if ctrl.swANorm && ~ctrl.swNegCut
%     % Modify derivatives of alpha for applying normalization
%     dAlpha = compose(repmat(-sum(dAlpha .* alpha,1) / sigma.anorm, ...
%         [size(alpha,1),1]),dAlpha,alpha);
% end
% 
% % calculate final derivatives of alpha with smoothness prior
% if ctrl.swPatOpt
%     % reshape alpha for convenience
%     alpha = reshape(alpha,[res,size(beta,2)]);
%     % intermediate results
%     drow = diff(alpha,1,1);
%     dcol = diff(alpha,1,2);
%     % calculate original smooth derivatives
%     dSmooth = - (diff(padarray(drow,[1,0]),1,1) + diff(padarray(dcol,[0,1]),1,2));
% %     % modify derivative for alpha-normalization
% %     if ctrl.swANorm && ~ctrl.swNegCut
% %         dSmooth = compose(repmat(-(sum(sum(drow.^2,1),2) ...
% %             + sum(sum(dcol.^2,1),2)) / sigma.anorm,[res,1]),dSmooth,alpha);
% %     end
%     % assemble final derivatives of alpha
%     dAlpha = dAlpha / sigma.noise^2 + reshape(dSmooth,size(dAlpha)) ...
%         / (2 * numel(alpha) * sigma.smpat^2 / numel(delta));
% end
% 
% % calculate final derivatives of phi with smoothness prior
% if ctrl.swTransOpt
%     phi  = reshape(phi,[res,size(theta,3)]);
%     dPhi = dPhi / sigma.noise^2 - reshape( ...
%         diff(wrapToPi(diff(padarray(phi,[1,0,0],'replicate'),1,1)),1,1) ...
%         + diff(wrapToPi(diff(padarray(phi,[0,1,0],'replicate'),1,2)),1,2),size(dPhi)) ...
%         / (2 * numel(phi) * sigma.smtrans^2 / numel(delta));
% end

% Normalize Gradients
normFactor = sqrt(sum(dAlpha(:).^2) + sum(dPhi(:).^2));
if ctrl.swPatOpt,  dAlpha = dAlpha / normFactor; end
if ctrl.swTransOpt, dPhi   = dPhi   / normFactor; end
end

%==========================================================================

% function S = sign(X)
% % SIGN function here is different from the system buildin one, in this
% % version when x equals to 0, 1 instead of 0 returns. This is required by
% % optimization procedure here.
% S = ones(size(X));
% S(X<0) = -1;
% end

function D = compose(S,U,alpha)
% COMPOSE function compose derivatives for alpha by stable and unstable
% parts. This function is created to deal with the discontinuity of
% abosolute value function in alpha-normalization.
assert(all(size(S)==size(U)),'two input matrix have to have same size!');
% initialize derivatives with value of negative region
D = U - S;
% modify derivatives in positive region
map = alpha > 0; D(map) = S(map) + U(map);
% modify derivatives in zero region
% - calculate index maps for decision tree
map0 = alpha == 0;
mapA = (S + U) <= 0;
mapB = U <= 0;
mapC = S >= U;
% - compose derivatives
map =  mapA & mapB & map0; D(map) = S(map) + U(map);
map = ~mapA & mapC & map0; D(map) = 0;
end

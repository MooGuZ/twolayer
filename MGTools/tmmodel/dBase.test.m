function [dalpha, dphi, normFactor] = ...
    dBase(alpha, phi, beta, theta, bias, delta, sigma, ctrl, res)
% TMDBASE calculate derivatives of basis functions in Transformation-Mask
% model.
%
% USAGE : [dalpha,dphi] = tmdbase(alpha,phi,beta,theta,bia,delta,sigma,res,option)
%
% MooGu Z. <hzhu@case.edu>
% June 17, 2015 - Version 0.00 : initial commit

% short path for easy case
if ~(ctrl.swPatOpt || ctrl.swTransOpt)
	dalpha = 0;
	dphi = 0;
	normFactor = 1;
	return;
end

% initialize gradients
if isa(alpha, 'gpuArray')
    dalpha = gpuArray.zeros(size(alpha));
    dphi   = gpuArray.zeros(size(phi));
else
    dalpha = zeros(size(alpha));
    dphi   = zeros(size(phi));
end

% constant factors in probability description
fnoise   = sigma.noise^2   * numel(delta);
fsmmask  = 2 * sigma.smpat^2  * numel(alpha);
fsmtrans = 2 * sigma.smtrans^2 * numel(phi);

% calculate NOISE part of gradient
for f = 1 : size(delta, 4)
    % <phase> as common term in gradient calculation of <alpha> and <phi>
    phase = bsxfun(@minus, phi, theta(:,:,:,f));
    
    % gradient of <alpha> in NOISE (likelihood)
    dalpha = dalpha - bsxfun(@times, delta(:,:,:,f), bsxfun(@plus, bias(:,:,:,f), ...
        sum(bsxfun(@times, beta(:,:,:,f), cos(phase)), 3)));
    
    % gradient of <phi> in NOISE (likelihood)
    dphi = dphi + bsxfun(@times, delta(:,:,:,f), sum(bsxfun(@times, alpha, ...
        beta(:,:,:,f)) .* sin(phase), 2));
end

% modify gradient of <alpha> in NOISE by introducing normalization
% - this operation need <all(alpha >= 0)> and need post process cut-off 
% - negative <alpha> to 0
if ctrl.swANorm && ctrl.swNegCut
	dalpha = bsxfun(@plus, dalpha, ...
    	sum(dalpha .* alpha, 1) / ctrl.anorm);
end	

% finalize gradient in NOISE
dalpha = dalpha / fnoise;
dphi   = dphi   / fnoise;

% reshape bases for computational convenience in SMOOTH
alpha = reshape(alpha, [res, size(beta, 2)]);
phi   = reshape(phi, [res, size(theta, 3)]);

% intermediate result
drow = diff(alpha, 1, 1);
dcol = diff(alpha, 1, 2);

% gradient of <alpha>
% - this operation need <all(alpha >= 0)> and need post process cut-off 
% - negative <alpha> to 0
dsmooth = -(diff(padarray(drow, [1,0]), 1, 1) ...
	+ diff(padarray(dcol, [0,1]), 1, 2));
	
if ctrl.swANorm && ctrl.swNegCut
	dsmooth = bsxfun(@minus, dsmooth, (sum(sum(drow.^2, 1), 2) ...
		+ sum(sum(dcol.^2, 1), 2)) / ctrl.anorm);
end

dalpha = dalpha + reshape(dsmooth / fsmmask, [prod(res), size(beta,2)]);

% gradient of <phi>
dphi = dphi - reshape( ...
    diff(padarray(wrapToPi(diff(phi, 1, 1)), [1,0]), 1, 1) + ...
    diff(padarray(wrapToPi(diff(phi, 1, 2)), [0,1]), 1, 2), ...
    prod(res), 1, size(theta, 3)) / fsmtrans;
	
% Normalization
normFactor = sqrt(sum(dalpha(:).^2) + sum(dphi(:).^2));
dalpha = dalpha / normFactor;
dphi   = dphi / normFactor;

if ~ctrl.swPatOpt, dalpha = 0; end
if ~ctrl.swTransOpt, dphi = 0; end

end

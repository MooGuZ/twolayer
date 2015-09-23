function [dbeta, dtheta, dbia, normFactor] = ...
    dCoefficient(alpha, phi, beta, theta, bia, delta, sigma, ffindex)
% TMDCOEF calculate derivatives of coefficients of the data in
% representation of Transformation-Mask model.
%
% USAGE : [dbeta,dtheta,dbia] = tmdcoef(alpha,phi,beta,theta,bia,delta,sigma)
%
% MooGu Z. <hzhu@case.edu>
% June 17, 2015 - Version 0.00 : initial commit

% initialize derivatives
if isa(alpha, 'gpuArray')
    dbia   = gpuArray.zeros(size(bia));
    dbeta  = gpuArray.zeros(size(beta));
    dtheta = gpuArray.zeros(size(theta));
else
    dbia   = zeros(size(bia));
    dbeta  = zeros(size(beta));
    dtheta = zeros(size(theta));
end

% difference of adjacent frame of <theta>
segDiff = wrapToPi(diff(theta, 1, 4));
segDiff(:, :, :, ffindex(2:end)-1) = 0;
segDiff = padarray(segDiff, [0,0,0,1]);

% constant factors of each probability term (noise, sparse, and slow)
fnoise  = sigma.noise^2 * numel(delta);
fsparse = numel(beta);
fslow   = sigma.slow^2 * numel(theta);

% calculate derivative frame by frame
for f = 1 : size(beta,4)
    % Common term in derivative of beta, theta and bia
    phase    = bsxfun(@minus,phi,theta(:,:,:,f));
    errorMap = bsxfun(@times,delta(:,:,:,f),alpha);
    % Derivatives
    dbia(:,:,:,f)   = -sum(errorMap,1) / fnoise;
    dbeta(:,:,:,f)  = -sum(bsxfun(@times,errorMap,cos(phase)),1) / fnoise ...
        + beta(:,:,:,f) ./ (beta(:,:,:,f).^2 + sigma.sparse^2) / fsparse;
    dtheta(:,:,:,f) = (segDiff(:,:,:,f+1) - segDiff(:,:,:,f)) / fslow ...
        - (beta(:,:,:,f) .* sum(bsxfun(@times,errorMap,sin(phase)),1)) / fnoise;
end

% Normalize Gradient
normFactor = sqrt(sum(dbia(:).^2) + sum(dbeta(:).^2) + sum(dtheta(:).^2));
dbia   = dbia   / normFactor;
dbeta  = dbeta  / normFactor;
dtheta = dtheta / normFactor;

end

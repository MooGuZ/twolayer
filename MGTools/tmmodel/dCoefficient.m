function [dBeta,dTheta,dBia,normFactor] = ...
    dCoefficient(alpha,phi,beta,theta,bia,delta,sigma,ffindex)
% initialize derivatives
dBia   = zeros(size(bia));
dBeta  = zeros(size(beta));
dTheta = zeros(size(theta));
if gpuDeviceCount ~= 0
    dBia   = gsingle(dBia);
    dBeta  = gsingle(dBeta);
    dTheta = gsingle(dTheta);
end
% Intermediate result of slow prior for theta
if gpuDeviceCount == 0
    segDiff = wrapToPi(diff(theta,1,4));
    segDiff(:,:,:,ffindex(2:end)-1) = 0;
    dSlow = -diff(padarray(segDiff,[0,0,0,1]),1,4);
else
    theta = reshape(theta,size(theta,2),size(theta,3),size(theta,4));
    segDiff = wrapToPi(diff(theta,1,3));
    segDiff(:,:,ffindex(2:end)-1) = 0;
    dSlow = reshape(-diff(padarray(segDiff,[0,0,1]),1,3), ...
                    [1,size(theta)]);
    theta = reshape(theta,[1,size(theta)]);
end
% Scale ratio for derivative of beta and theta
ratio = numel(delta) / numel(beta);
% calculate derivative frame by frame
for f = 1 : size(beta,4)
    % Common term in derivative of beta, theta and bia
    phase    = bsxfun(@minus,phi,theta(:,:,:,f));
    errorMap = bsxfun(@times,delta(:,:,:,f),alpha);
    % Derivatives
    dBia(:,:,:,f)   = -sum(errorMap,1) / sigma.noise^2;
    dBeta(:,:,:,f)  = -sum(bsxfun(@times,errorMap,cos(phase)),1) / sigma.noise^2 + ...
        (ratio * beta(:,:,:,f)) ./ (beta(:,:,:,f).^2 + sigma.sparse^2);
    dTheta(:,:,:,f) = -(beta(:,:,:,f) .* sum(bsxfun(@times,errorMap, ...
        sin(phase)),1)) / sigma.noise^2 + ratio * dSlow(:,:,:,f) / sigma.slow^2;
end
% Normalize Gradient
normFactor = sqrt(sum(dBia(:).^2) + sum(dBeta(:).^2) + sum(dTheta(:).^2));
dBia   = dBia   / normFactor;
dBeta  = dBeta  / normFactor;
dTheta = dTheta / normFactor;
end

function obj = objFunc(alpha,phi,beta,theta,~,delta,sigma,ffindex,res)
npixel   = size(alpha,1);
npattern = size(beta,2);
ntrans   = size(phi,3);
nframe   = size(theta,4);
nrow     = res(1);
ncol     = res(2);
% Difference of each segment along time axis of theta
% if gpuDeviceCount == 0
    segDiff = wrapToPi(diff(theta,1,4));
    segDiff(:,:,:,ffindex(2:end)-1) = 0;
% else
%     theta = reshape(theta,npattern,ntrans,nframe);
%     segDiff = wrapToPi(diff(theta,1,3));
%     segDiff(:,:,ffindex(2:end)-1) = 0;
%     segDiff = reshape(segDiff,1,npattern,ntrans,nframe-1);
% end
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
obj.value   = obj.noise + obj.sparse + obj.slow + obj.smpat + obj.smtrans;
end

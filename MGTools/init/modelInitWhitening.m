function [m,p] = modelInitWhitening(m,p)
% MODELINITWHITENING calculate parameters that needed by whitening process
%
%   [m,p] = modelInitWhitening(m,p)
%
% MooGu Z. <hzhu@case.edu>
% June 17, 2014 - Version 0.1

% ======== Load Data =======
Data = zeros(m.patch_sz^2,p.data.scope*p.data.nframe,'double');
if p.data.scope >= p.data.quantity
    for i = 1 : p.data.scope
        Data(:,(i-1)*p.data.nframe+1:i*p.data.nframe) ...
            = gif2anim([p.data.path,p.data.nameList{i}],1:p.data.nframe);
    end
else  
    % Random Number Genereating
    s = randi(p.data.quantity,p.data.scope,1);
    % Load GIF Files
    for i = 1 : p.data.scope
        Data(:,(i-1)*p.data.nframe+1:i*p.data.nframe) ...
            = gif2anim([p.data.path,p.data.nameList{s(i)}],1:p.data.nframe);
    end
end

% ======= Parameter Estimation =======
m.imageMean = mean(Data,2);
Data = bsxfun(@minus,Data,m.imageMean); % Remove Bias
m.pixel_variance = var(Data(:));
m.pixel_noise_variance = p.whitening.pixel_noise_fractional_variance * ...
    m.pixel_variance;
% ======= Eigen Value Analysis =======
C = Data * Data' / (p.data.scope*p.data.nframe); clear Data
[eigVec,eigVal] = eig(C);
[eigVal,index]  = sort(diag(eigVal),'descend');
eigVec = eigVec(:,index);
% Save in Model
m.imageEigVals = sparse(diag(eigVal));
m.imageEigVecs = eigVec;
% ======= Determine Dimension of Whitening Signal =======
varCutoff = p.whitening.pixel_noise_variance_cutoff_ratio * m.pixel_noise_variance;
m.M = sum(eigVal > varCutoff);

eigVal = eigVal(1:m.M);
eigVec = eigVec(:,1:m.M);
D = diag(real(eigVal.^(-0.5)));

iRolloff = sum(eigVal > varCutoff * p.whitening.X_noise_fraction);
m.I_noise_factors = ones(m.M,1);
m.I_noise_factors(iRolloff+1 : end) = .5*(1+cos(linspace(0,pi,m.M-iRolloff)));
m.I_noise_factors = m.I_noise_factors / p.whitening.X_noise_var;
m.I_noise_vars = 1 ./ m.I_noise_factors;

% whitening transform
m.whitenMatrix = D * eigVec';
m.dewhitenMatrix = eigVec / D;
m.zerophaseMatrix = eigVec * D * eigVec';

end
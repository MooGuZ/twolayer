function Model = tpmodel(videoPath,varargin)
% TPMODEL learn Transformation and Pattern separated base function sets
%
% MODEL = TPMODEL(VIDEOPATH,SETTINGS...) learn transformation and pattern
% bases sets from videos from VIDEOPATH with specific SETTINGS
%
% Available Settings:
% -----------------------------------------------------------------------
% NAME                    |DESCRIPTION                           |DEFAULT
% ------------------------|--------------------------------------|-------
% Model                   |previous model                        | []
% nPattern                |number of pattern bases               | 13
% nTrans                  |number of transformation bases        |  2
% nEpoch                  |number of optimization rounds         | 13
% nAdapt                  |number of base adapting steps         | 30
% nInfer                  |number of coefficient infering steps  | 30
% InitialStepRatio        |initial step in search space ratio    | 0.01
% Accuracy                |optimization precision requirement    | 0.0001
% NoisePrior              |sigma of noise distribution           | 0.1
% SparsePrior             |sigma of sparseness distribution      | 1
% SlowPrior               |sigma of slowness distribution        | 2 PI
% PatternBaseSmoothPrior  |sigma of smoothness of pattern base   | 1
% TransformBaseSmoothPrior|sigma of smoothness of transform base | 2 PI
% OptimizePatternBase     |swicher of pattern base optimization  | TRUE
% OptimizeTransformBase   |switcher of transformbase optimization| TRUE
% -----------------------------------------------------------------------
%
% MooGu Z. <hzhu@case.edu>
% Jan 21, 2015 - Version 0.1
% Feb 07, 2015 - Version 0.2

% Input Parser
% ------------
% initialize input parser
p = inputParser;
% define input parameters
p.addRequired('videoPath', @(x) exist(x,'file'));
p.addParamValue('Model',[], @(x) ...
    isstruct(x) && isfield(x,'transBase') && isfield(x,'patBase') ...
    && isfield(x,'transCoefficient') && isfield(x,'patCoefficient') ...
    && isfield(x,'bia'));
p.addParamValue('nPattern', 13, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nTrans', 2, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nEpoch', 13, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nAdapt', 30, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nInfer', 30, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('InitialStepRatio', 1e-2, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('Accuracy', 1e-4, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('NoisePrior',1e-1, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('SparsePrior', 1e0, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('SlowPrior', 2e0 * pi, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('PatternBaseSmoothPrior', 1e0, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('TransformBaseSmoothPrior', 2e0 * pi, @(x) ...
    isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('OptimizePatternBase', true, @(x) ...
    islogical(x) && isscalar(x));
p.addParamValue('OptimizeTransformBase', true, @(x) ...
    islogical(x) && isscalar(x));
% parse input
p.parse(videoPath,varargin{:});
% initialize parameter according to input argument parsing result
% 1. workflow parameters
npattern = p.Results.nPattern;
ntrans   = p.Results.nTrans;
nEpoch   = p.Results.nEpoch;
nAdapt   = p.Results.nAdapt;
nInfer   = p.Results.nInfer;
stepInit = p.Results.InitialStepRatio;
accuracy = p.Results.Accuracy;
% 2. statistic priors
sigma.noise   = p.Results.NoisePrior;
sigma.sparse  = p.Results.SparsePrior;
sigma.slow    = p.Results.SlowPrior;
sigma.smpat   = p.Results.PatternBaseSmoothPrior;
sigma.smtrans = p.Results.TransformBaseSmoothPrior;
% 3. optimization switcher
swPatOpt   = p.Results.OptimizePatternBase;
swTransOpt = p.Results.OptimizeTransformBase;
% 4. previous model
Model = p.Results.Model;

% Check availability of GPU
swGPU = (gpuDeviceCount ~= 0);

% Generating File List
switch exist(videoPath,'file')
    case {2}
        flist = {videoPath};
    case {7}
        flist = dir([videoPath,'/*.gif']);
        flist = strcat([videoPath,'/'],{flist(:).name});
    otherwise
        error('Input argument #1 is not refering any folder or GIF file!');
end

% Generate Collection of Animations and list of first/last frame index
ffindex = 1; % First Frame Index
[v,animRes] = gif2anim(flist{1});
for i = 2 : numel(flist)
    ffindex = [ffindex,size(v,2)+1];
    v = [v,gif2anim(flist{i})];
end

% GPU Enabling
if swGPU
	v = gsingle(v);
end

% Quantities of Data
[npixel,nframe] = size(v);

% Initialize Base Funtion Set and Parameters
% In this program all the data arranged in the coordinate system with axes
% <PIXEL * FRAME * PATTERN * TRANSFORMATION>
if isempty(Model)
    phi    = wrapToPi(pi * randn(npixel,1,1,ntrans));
    alpha  = rand(npixel,1,npattern);
    theta  = wrapToPi(pi * randn(1,nframe,npattern,ntrans));
    beta   = randn(1,nframe,npattern,ntrans);
    bia    = rand(1,nframe,npattern,1);
else
    [alpha,phi,beta,theta,bia] = tpmodel2data(Model);
    [~,npattern] = size(alpha);
    [~,ntrans]   = size(phi);
    assert(npixel==size(alpha,1), ...
        'Model must have the same pixel number as animation.');
    assert(nframe==size(beta,2), ...
        'Model must have the same frame number as animation.');
end

% GPU Enabling
if swGPU
	phi = gsingle(phi);
	alpha = gsingle(alpha);
	theta = gsingle(theta);
	beta = gsingle(beta);
	bia = gsingle(bia);
end	

% Initialize Objective Value Records
objRec = zeros(7,2*nEpoch+1);

% Calculate initial step size for inference and adaption
stepInitInfer = stepInit * sqrt(npattern*nframe*(1+8*pi^2*ntrans));
stepInitAdapt = 0;
if swPatOpt
    stepInitAdapt = stepInitAdapt + stepInit * sqrt(npixel*npattern); 
end
if swTransOpt
    stepInitAdapt = stepInitAdapt + stepInit * sqrt(npixel*4*pi^2*ntrans);
end

% E-M Algo
delta = v - genModel(alpha,phi,beta,theta,bia);
objective = objFunc(alpha,phi,beta,theta,bia,delta,sigma,ffindex,animRes);
for epoch = 1 : nEpoch
    % Infering Optimal Theta
    [beta,theta,bia,delta,objective,niter] = ...
        inferTPM(nInfer,alpha,phi,beta,theta,bia,delta,objective, ...
            v,sigma,ffindex,animRes,stepInitInfer,accuracy);
    % Show information
    disp(['Objective Value after infering process of EPOCH[', ...
        num2str(epoch),'] >> ',num2str(objective.value), ...
        ' (',num2str(niter),' cycles)']);
    % Records Objective Values
    objRec(:,2*epoch-1) = [niter,objective.value,objective.noise, ...
        objective.sparse,objective.slow,objective.smpat,objective.smtrans]';
    
    % Adapting Complex Base Function
    [alpha,phi,delta,objective,niter] = ...
        adaptTPM(nAdapt,alpha,phi,beta,theta,bia,delta,objective, ...
            v,sigma,ffindex,animRes,stepInitAdapt,accuracy, ...
            swPatOpt,swTransOpt);
    % Show information
    disp(['Objective Value after adapting process of EPOCH[', ...
        num2str(epoch),'] >> ',num2str(objective.value), ...
        ' (',num2str(niter),' cycles)']);
    % Records Objective Values
    objRec(:,2*epoch) = [niter,objective.value,objective.noise, ...
        objective.sparse,objective.slow,objective.smpat,objective.smtrans]';
    
    % Normalize alpha
    [alpha,phi,beta,theta,bia,delta,objective] = ...
        normalizeAlpha(alpha,phi,beta,theta,bia,delta,objective, ...
            animRes,v,sigma,ffindex);
end
% For case nEpoch == 0
if isempty(epoch), epoch = 0; end

% Inference for the final bases
[beta,theta,bia,~,objective,niter] = ...
    inferTPM(nInfer,alpha,phi,beta,theta,bia,delta,objective, ...
        v,sigma,ffindex,animRes,stepInitInfer,accuracy);
% Show information
disp(['Objective Value after final infering process >> ', ...
    num2str(objective.value),' (',num2str(niter),' cycles)']);
% Records Objective Values
objRec(:,2*epoch+1) = [niter,objective.value,objective.noise, ...
        objective.sparse,objective.slow,objective.smpat,objective.smtrans]';

% Data of Model
Model.transBase        = reshape(phi,npixel,ntrans);
Model.patBase          = reshape(alpha,npixel,npattern);
Model.transCoefficient = permute(theta,[3,4,2,1]);
Model.patCoefficient   = permute(beta,[3,4,2,1]);
Model.bia              = permute(bia,[3,2,1,4]);
% Information of Model
Model.sigma            = sigma;
if isfield(Model,'objRec')
    Model.objRec = [Model.objRec,objRec];
else
    Model.objRec = objRec;
end
% Recovered Animation
Model.recoveredAnim = genModel(alpha,phi,beta,theta,bia);
% Save Index of first frame
Model.ffindex = ffindex;

% Tranform data from GPU format to CPU format
if swGPU
	Model.transBase = double(Model.transBase);
	Model.patBase   = double(Model.patBase);
	Model.transCoefficient = double(Model.transCoefficient);
	Model.patCoefficient   = double(Model.patCoefficient);
	Model.bia = double(Model.bia);
end

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

function [beta,theta,bia,delta,obj,i] = ...
    inferTPM(niter,alpha,phi,beta,theta,bia,delta,obj, ...
        v,sigma,ffindex,resolution,stepInit,stepMin)
if niter < 1 
    i = 0; return
end

step = stepInit;
for i = 1 : niter
    [dBeta,dTheta,dBia] = ...
        dCoefficient(alpha,phi,beta,theta,bia,delta,sigma,ffindex);
    newBeta  = beta - step * dBeta;
    newTheta = wrapToPi(theta - step * dTheta);
    newBia   = bia - step * dBia;
    newDelta = v - genModel(alpha,phi,newBeta,newTheta,newBia);
    newObj   = ...
        objFunc(alpha,phi,newBeta,newTheta,newBia,newDelta, ...
            sigma,ffindex,resolution);
    while(newObj.value > obj.value)
        step = step / 2;
        if step < stepMin, break; end
        newBeta  = beta - step * dBeta;
        newTheta = wrapToPi(theta - step * dTheta);
        newBia   = bia - step * dBia;
        newDelta = v - genModel(alpha,phi,newBeta,newTheta,newBia);
        newObj   = ...
            objFunc(alpha,phi,newBeta,newTheta,newBia,newDelta, ...
                sigma,ffindex,resolution);
    end
    if step < stepMin, break; end
    beta  = newBeta;
    theta = newTheta;
    bia   = newBia;
    delta = newDelta;
    obj   = newObj;
end
end

function [alpha,phi,delta,obj,i] = ...
    adaptTPM(niter,alpha,phi,beta,theta,bia,delta,obj, ...
        v,sigma,ffindex,resolution,stepInit,stepMin, ...
        swAlpha,swPhi)
if niter < 1 || ~(swAlpha || swPhi)
    i = 0; return
end

step = stepInit;
for i = 1 : niter
    [dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia,delta, ...
        sigma,resolution,swAlpha,swPhi);
    if swAlpha
        newAlpha = alpha - step * dAlpha; 
    else
        newAlpha = alpha; 
    end
    if swPhi
        newPhi = wrapToPi(phi - step * dPhi);
    else
        newPhi = phi;
    end
    newDelta = v - genModel(newAlpha,newPhi,beta,theta,bia);
    newObj   = ...
        objFunc(newAlpha,newPhi,beta,theta,bia,newDelta, ...
            sigma,ffindex,resolution);
    while(newObj.value > obj.value)
        step = step / 2;
        if step < stepMin, break; end
        if swAlpha
            newAlpha = alpha - step * dAlpha;
        else
            newAlpha = alpha;
        end
        if swPhi
            newPhi = wrapToPi(phi - step * dPhi);
        else
            newPhi = phi;
        end
        newDelta = v - genModel(newAlpha,newPhi,beta,theta,bia);
        newObj   = ...
            objFunc(newAlpha,newPhi,beta,theta,bia,newDelta, ...
                sigma,ffindex,resolution);
    end
    if step < stepMin, break; end
    alpha = newAlpha;
    phi   = newPhi;
    delta = newDelta;
    obj   = newObj;
end
end

function [alpha,phi,beta,theta,bia,delta,obj] = ...
    normalizeAlpha(alpha,phi,beta,theta,bia,delta,obj,res,v,sigma,ffindex)
zeroDecisionBound = 1e-2;
% Scale factors of alpha
a = max(abs(alpha));
% Positive map of alpha
P = (alpha < 0);
% Set X and I
setX = any(P,3);
setI = any(P,1);
% Normalization and Transformation by Scale Factor
if any(abs(a-1) > zeroDecisionBound)
    alpha = bsxfun(@rdivide,alpha,a);
    bia   = bsxfun(@times,bia,a);
    beta  = bsxfun(@times,beta,a);
end
% Normalization and Transformation by Positive Restraint
if any(setX) || any(setI)
    alpha = abs(alpha);
%     % Approximate Compensation of Absolute Operator
%     if sum(setX)*sum(setI) > sum(~setX)*sum(~setI)
%         phi(setX(:),:,:,:) = wrapToPi(phi(setX(:),:,:,:) + pi);
%         theta(:,:,~setI(:),:) = wrapToPi(theta(:,:,~setI(:),:) + pi);
%     end
%     % Send warning for disability of compansation normalization
%     if (sum(~setX) * sum(~setI) ~= 0) || ...
%             any(all([permute(any(abs(bia)>zeroDecisionBound,2),[3,1,2,4]), ...
%             setI(:)],2))
%         disp('[Warning] The effects of alpha normalization cannot be fully compansated!');
%     end
    % Recalculate Objective Value
    delta = v - genModel(alpha,phi,beta,theta,bia);
    obj   = objFunc(alpha,phi,beta,theta,bia,delta,sigma,ffindex,res);
    disp(['Objective Value after normalization of alpha >> ', ...
        num2str(obj.value)]);
end

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

function [dBeta,dTheta,dBia] = ...
    dCoefficient(alpha,phi,beta,theta,~,delta,sigma,ffindex)
% Intermediate result of slow prior for theta
segDiff = wrapToPi(diff(theta,1,2));
segDiff(:,ffindex(2:end)-1,:,:) = 0;
% Common term in derivative of beta, theta and bia
phase    = bsxfun(@minus,phi,theta);
errorMap = bsxfun(@times,delta,alpha);
% Scale ration for derivative of beta and theta
ratio    = numel(delta) / numel(beta);
% Derivatives
dBia   = -sum(errorMap,1) / sigma.noise^2;
dBeta  = -sum(bsxfun(@times,errorMap,cos(phase)),1) / sigma.noise^2 + ...
    (ratio * beta) ./ (beta.^2 + sigma.sparse^2);
dTheta = -(beta .* sum(bsxfun(@times,errorMap,sin(phase)),1)) / sigma.noise^2 + ...
    ratio * (-diff(padarray(segDiff,[0,1,0,0]),1,2)) / sigma.slow^2;
% Normalize Gradient
normFactor = sqrt(sum(dBia(:).^2) + sum(dBeta(:).^2) + sum(dTheta(:).^2));
dBia   = dBia   / normFactor;
dBeta  = dBeta  / normFactor;
dTheta = dTheta / normFactor;
end

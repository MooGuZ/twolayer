function [model,rec,video] = tmmodel(video,varargin)
% TPMODEL learn Transformation and Pattern separated base function sets
%
% MODEL = TPMODEL(VIDEO,SETTINGS...) learn transformation and pattern
% bases sets from videos from VIDEO with specific SETTINGS, VIDEO can be a
% path of file or folder, it also can be a structure contraining
% corresponding information.
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
% NegativeCutProbability  |probability of applying negative cut  | 0.7
% OptimizePatternBase     |swicher of pattern base optimization  | TRUE
% OptimizeTransformBase   |switcher of transformbase optimization| TRUE
% AlphaNormalization      |switcher of normalizing alpha         | FALSE
% Verbose                 |verbose level of shown information    | 2
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
p.addRequired('video', @(x) ...
    isstruct(x) || exist(x,'file'));
p.addParamValue('model',[], @(x) ...
    isstruct(x) && isfield(x,'phi') && isfield(x,'alpha') ...
    && isfield(x,'theta') && isfield(x,'beta') && isfield(x,'bia'));
p.addParamValue('nPattern', 13, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nTrans', 2, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nEpoch', 13, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nAdapt', 30, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nInfer', 30, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('InitialStepRatio', 1e-2, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('Accuracy', 1e-4, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('NoisePrior',1e-1, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('SparsePrior', 1e0, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('SlowPrior', 2e0 * pi, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('PatternBaseSmoothPrior', 1e0, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('TransformBaseSmoothPrior', 2e0 * pi, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('NegativeCutProbability', 0.7, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
p.addParamValue('OptimizePatternBase', true, ...
    @(x) islogical(x) && isscalar(x));
p.addParamValue('OptimizeTransformBase', true, ...
    @(x) islogical(x) && isscalar(x));
p.addParamValue('AlphaNormalization', false, ...
    @(x) islogical(x) && isscalar(x));
p.addParamValue('Verbose', 2, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x));
% parse input
p.parse(video,varargin{:});
% initialize parameter according to input argument parsing result
% 1. workflow parameters
npattern = p.Results.nPattern;
ntrans   = p.Results.nTrans;
nEpoch   = p.Results.nEpoch;
nAdapt   = p.Results.nAdapt;
nInfer   = p.Results.nInfer;
stepInit = p.Results.InitialStepRatio;
% 2. statistic priors
sigma.noise   = p.Results.NoisePrior;
sigma.sparse  = p.Results.SparsePrior;
sigma.slow    = p.Results.SlowPrior;
sigma.smpat   = p.Results.PatternBaseSmoothPrior;
sigma.smtrans = p.Results.TransformBaseSmoothPrior;
% 3. probability of applying negative-cut
ctrl.probNegCut = p.Results.NegativeCutProbability;
% 4. optimization switcher
ctrl.swPatOpt   = p.Results.OptimizePatternBase;
ctrl.swTransOpt = p.Results.OptimizeTransformBase;
% 5. alpha normalization switcher
ctrl.swANorm  = p.Results.AlphaNormalization;
% 6. optimization accuracy
ctrl.accuracy = p.Results.Accuracy;
% 7. verbose level of output information
ctrl.verbose  = p.Results.Verbose;
% 8. previous model
model = p.Results.model;

% Check availability of GPU
% swGPU = (gpuDeviceCount ~= 0);
swGPU = false;

% prepare video data
if isstruct(video)
    % create alias
    v = video.v;
    animRes = video.res;
    ffindex = video.ffindex;
else
    % read animation files
    [v,ffindex,animRes] = dataPrepare(video);
end

% Quantities of Data
[npixel,nframe] = size(v);

% Initialize Base Funtion Set and Parameters
% In this program all the data arranged in the coordinate system with axes
% <PIXEL * FRAME * PATTERN * TRANSFORMATION>
if isempty(model)
    phi   = wrapToPi(pi * randn(npixel,1,ntrans,1));
    alpha = rand(npixel,npattern,1,1);
    theta = wrapToPi(pi * randn(1,npattern,ntrans,nframe));
    beta  = randn(1,npattern,ntrans,nframe);
    bia   = rand(1,npattern,1,nframe);
    % write probabilistic parameter to model
    model.sigma = sigma;
else
    % renew quantity information
    npattern = size(model.alpha,2);
    ntrans   = size(model.phi,2);
    % reshape model parameters for optimization
    [alpha,phi,beta,theta,bia] = m2p(model);
    % get probabilistic setting
    if isfield(model,'sigma')
        sigma = model.sigma;
        disp('Paramter SIGMA is set by input model structure!');
    end
    % get functional paramters
    if isfield(model,'ctrl')
        ctrl  = model.ctrl;
        disp('Parameter CTRL is set by input model structure!');
    end
end
% set up normalized value of alpha
if ctrl.swANorm && ~isfield(ctrl,'anorm')
    % calculate sum-value of each pattern base
    svalue = sum(abs(alpha),1);
    % calculate normalize value of alpha
    ctrl.anorm = mean(svalue,2);
    % calculate scaling ratio of each pattern base
    sratio = svalue / ctrl.anorm;
    % normalize alpha with ocpratio
    alpha  = bsxfun(@rdivide,alpha,sratio);
    % compansating beta and bias
    beta   = bsxfun(@times,beta,sratio);
    bia    = bsxfun(@times,bia,sratio);
end
% reshape video for calculating convenience
v = reshape(v,[npixel,1,1,nframe]);

% GPU Enabling
if swGPU
    phi   = gsingle(phi);
    alpha = gsingle(alpha);
    theta = gsingle(theta);
    beta  = gsingle(beta);
    bia   = gsingle(bia);
    v     = gsingle(v);
end	

% Initialize Objective Value Records
objRec.n = zeros(1,2*nEpoch+1);
objRec.v = repmat(struct('noise',0,'sparse',0,'slow',0, ...
    'smpat',0,'smtrans',0,'value',0),1,2*nEpoch+1);

% Calculate initial step size for inference and adaption
ctrl.inferInitStep = stepInit * sqrt(npattern*nframe*(1+8*pi^2*ntrans));
ctrl.adaptInitStep = 0;
if ctrl.swPatOpt
    ctrl.adaptInitStep = ctrl.adaptInitStep + stepInit * sqrt(npixel*npattern); 
end
if ctrl.swTransOpt
    ctrl.adaptInitStep = ctrl.adaptInitStep + stepInit * sqrt(npixel*4*pi^2*ntrans);
end

% E-M Algo
delta = v - genmodel(alpha,phi,beta,theta,bia);
objective = objFunc(alpha,phi,beta,theta,bia,delta,sigma,ffindex,animRes);
for epoch = 1 : nEpoch
    % Infering Optimal Theta
    [beta,theta,bia,delta,objective,niter,ctrl] = ...
        inferGD(nInfer,alpha,phi,beta,theta,bia,delta,objective, ...
        sigma,ctrl,v,ffindex,animRes);
    % Records Objective Values
    objRec.n(2*epoch-1) = niter;
    objRec.v(2*epoch-1) = objective;
    % Show information
    if ctrl.verbose >= 2
        disp(['Objective Value after infering process of EPOCH[', ...
              num2str(epoch),'] >> ',num2str(objective.value), ...
              ' (',num2str(niter),' cycles)']);
    end
    
    % Adapting Complex Base Function
    [alpha,phi,delta,objective,niter,ctrl] = ...
        adaptGD(nAdapt,alpha,phi,beta,theta,bia,delta,objective, ...
            sigma,ctrl,v,ffindex,animRes);
    % Records Objective Values
    objRec.n(2*epoch) = niter;
    objRec.v(2*epoch) = objective;
    % Show information
    if ctrl.verbose >= 2
        disp(['Objective Value after adapting process of EPOCH[', ...
              num2str(epoch),'] >> ',num2str(objective.value), ...
              ' (',num2str(niter),' cycles)']);
    end
    
    % Normalize alpha
    [alpha,phi,beta,theta,bia,delta,objective] = ...
        normalizeAlpha(alpha,phi,beta,theta,bia,delta,objective, ...
            sigma,ctrl,v,ffindex,animRes,ctrl.verbose);
end
% For case nEpoch == 0
if isempty(epoch), epoch = 0; end

% Inference for the final bases
[beta,theta,bia,~,objective,niter,ctrl] = ...
    inferGD(nInfer,alpha,phi,beta,theta,bia,delta,objective, ...
        sigma,ctrl,v,ffindex,animRes);
% Records Objective Values
objRec.n(2*epoch+1) = niter;
objRec.v(2*epoch+1) = objective;
% Show information
if ctrl.verbose >= 1
    disp(['Objective Value after final infering process >> ', ...
          num2str(objective.value),' (',num2str(sum(objRec.n)),' cycles)']);
end

% Tranform data from GPU format to CPU format
if swGPU
    phi   = double(phi);
    alpha = double(alpha);
    theta = double(theta);
    beta  = double(beta);
    bia   = double(bia);
    v     = double(v);
end
    
% Data of Model
model = p2m(alpha,phi,beta,theta,bia,model);
% Information of Model
model.sigma = sigma;
model.ctrl  = ctrl;
model.obj   = objective;
if isfield(model,'objRec')
    model.objRec = objRecComb(model.objRec,objRec);
else
    for i = 2 : numel(objRec.n)
        objRec.n(i) = objRec.n(i) + objRec.n(i-1);
    end
    model.objRec = objRec;
end
% Data of Animation
if ~isstruct(video)
    clear('video')
    video.v = reshape(v,[npixel,nframe]);
    video.ffindex = ffindex;
    video.res = animRes;
end
% reconstruct animation by generative model
rec = reshape(genmodel(alpha,phi,beta,theta,bia),[npixel,nframe]);

end

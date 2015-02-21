function [model,video] = tpmodel(video,varargin)
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
p.addRequired('video', @(x) ...
    isstruct(x) || exist(x,'file'));
p.addParamValue('model',[], @(x) ...
    isstruct(x) && isfield(x,'phi') && isfield(x,'alpha') ...
    && isfield(x,'theta') && isfield(x,'beta') && isfield(x,'bia'));
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
p.parse(video,varargin{:});
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
model = p.Results.model;

% Check availability of GPU
swGPU = (gpuDeviceCount ~= 0);

% prepare video data
if isstruct(video)
    % create alias
    v = video.v;
    animRes = video.res;
    ffindex = video.ffindex;
else
    % Generating File List
    switch exist(video,'file')
        case {2}
            flist = {video};
        case {7}
            flist = dir([video,'/*.gif']);
            flist = strcat([video,'/'],{flist(:).name});
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
else
    % renew quantity information
    npattern = size(model.alpha,2);
    ntrans   = size(model.phi,2);
    % reshape model parameters for optimization
    [alpha,phi,beta,theta,bia] = m2p(model);
    % get probabilistic setting
    sigma = model.sigma;
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
delta = v - genmodel(alpha,phi,beta,theta,bia);
objective = objFunc(alpha,phi,beta,theta,bia,delta,sigma,ffindex,animRes);
for epoch = 1 : nEpoch
    % Infering Optimal Theta
    [beta,theta,bia,delta,objective,niter] = ...
        inferGD(nInfer,alpha,phi,beta,theta,bia,delta,objective, ...
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
        adaptGD(nAdapt,alpha,phi,beta,theta,bia,delta,objective, ...
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
    inferGD(nInfer,alpha,phi,beta,theta,bia,delta,objective, ...
        v,sigma,ffindex,animRes,stepInitInfer,accuracy);
% Show information
disp(['Objective Value after final infering process >> ', ...
    num2str(objective.value),' (',num2str(niter),' cycles)']);
% Records Objective Values
objRec(:,2*epoch+1) = [niter,objective.value,objective.noise, ...
        objective.sparse,objective.slow,objective.smpat,objective.smtrans]';

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
model.obj   = objective;
if isfield(model,'objRec')
    model.objRec = [model.objRec,objRec];
else
    model.objRec = objRec;
end
% Data of Animation
if ~isstruct(video)
    clear('video')
    video.v = reshape(v,[npixel,nframe]);
    video.ffindex = ffindex;
    video.res = animRes;
end
video.rec = reshape(genmodel(alpha,phi,beta,theta,bia),[npixel,nframe]);

end


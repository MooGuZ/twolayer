function model = tmminit(varargin)
% TMMINIT initialize a model for Transform-Mask model
%
% MODEL = TMMINIT(SETTINGS...) create a model with 7 mask bases and 7
% transform bases randomly generated. All other default values check
% following table.
%
% Available Settings:
% -----------------------------------------------------------------------
% NAME                    |DESCRIPTION                           |DEFAULT
% ------------------------|--------------------------------------|-------
% nPixel                  |number of pixels                      | 1024
% nMask                   |number of mask bases                  | 7
% nTrans                  |number of transformation bases        | 7
% MBase                   |specified mask bases                  | []
% TBase                   |specified transformation bases        | []
% Accuracy                |optimization precision requirement    | 0.0001
% NoisePrior              |sigma of noise distribution           | 1
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
% Apr 9, 2015

% Input Parser
% ------------
% initialize input parser
p = inputParser;
% define input parameters
p.addParamValue('nPixel', 1024, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nMask', 7, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('nTrans', 7, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && (floor(x) == x));
p.addParamValue('MBase', [], ...
    @(x) isnumeric(x) && isreal(x) && ismatrix(x));
p.addParamValue('TBase', [], ...
    @(x) isnumeric(x) && isreal(x) && ismatrix(x));
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
p.parse(varargin{:});
% initialize parameter according to input argument parsing result
% 1. workflow parameters
npixel = p.Results.nPixel;
nmask  = p.Results.nMask;
ntrans = p.Results.nTrans;
% 2. statistic priors
model.sigma.noise   = p.Results.NoisePrior;
model.sigma.sparse  = p.Results.SparsePrior;
model.sigma.slow    = p.Results.SlowPrior;
model.sigma.smpat   = p.Results.PatternBaseSmoothPrior;
model.sigma.smtrans = p.Results.TransformBaseSmoothPrior;
% 3. probability of applying negative-cut
model.ctrl.probNegCut = p.Results.NegativeCutProbability;
% 4. optimization switcher
model.ctrl.swPatOpt   = p.Results.OptimizePatternBase;
model.ctrl.swTransOpt = p.Results.OptimizeTransformBase;
% 5. alpha normalization switcher
model.ctrl.swANorm  = p.Results.AlphaNormalization;
% 6. optimization accuracy
model.ctrl.accuracy = p.Results.Accuracy;
% 7. verbose level of output information
model.ctrl.verbose  = p.Results.Verbose;
% 8. base functions
model.alpha = p.Results.MBase;
model.phi   = p.Results.TBase;

% randomly initialize mask bases
if isempty(model.alpha)
    model.alpha = rand(npixel,nmask,1,1);
end
% randomly initialize transform bases
if isempty(model.phi)
    model.phi = wrapToPi(pi * randn(npixel,1,ntrans,1));
end 

end
% This file aims at debuging the programs of transform-pattern model by
% creating sample bases and responds according to the distribution. This
% file is organized by cell units.
%
% MooGu Z. <hzhu@case.edu>
% Feb 19, 2015 - Version 0.1

clear; close all; rng shuffle; clc

%% fundamental settings
% number of bases
npattern = 2;
ntrans   = 1;
% size of frame
frmsz = [32,32];
% number of frams
nframe = 240;
% model noise
mnoise = 0.1;
% smoothness target value
smtgt = 0.1;
% parameters of prior and likelihood
sigma.noise   = 0.1;
sigma.sparse  = 1;
sigma.slow    = 2*pi;


%% generate reference bases

% initilialize bases
ref.alpha = zeros(prod(frmsz),npattern);
ref.phi   = zeros(prod(frmsz),ntrans);

ref.alpha(:,1) = reshape(pbaseGen([.5,.5],13,frmsz(1)),[prod(frmsz),1]);
ref.alpha(:,2) = reshape(pbaseGen([-.3,-.4],7,frmsz(1)),[prod(frmsz),1]);
ref.phi        = reshape(tbaseGen([1,1],3,frmsz(1)),[prod(frmsz),1]);

% estimate smoothness parameter
alpha = reshape(ref.alpha,[frmsz,npattern]);
phi   = reshape(ref.phi,[frmsz,ntrans]);
sigma.smpat = sqrt((sum(sum(sum(diff(alpha,1,1).^2))) ....
    + sum(sum(sum(diff(alpha,1,2).^2)))) ...
    / (2 * numel(alpha) * smtgt));
sigma.smtrans = sqrt((sum(sum(sum(wrapToPi(diff(phi,1,1)).^2))) ....
    + sum(sum(sum(wrapToPi(diff(phi,1,2)).^2)))) ...
    / (2 * numel(phi) * smtgt));

%% generate responds randomly

% save sigma setting to reference model
ref.sigma = sigma;

% beta, follow cauchy distribution
ref.beta = rand(npattern,ntrans,nframe);
ref.beta = sigma.sparse * tan(pi * (ref.beta - .5));

% theta, start from one frame and spread to all frames
% follow Gaussian distribution
ref.theta = sigma.slow * randn(npattern,ntrans,nframe);
for f = 2 : nframe
    ref.theta(:,:,f) = ref.theta(:,:,f-1) + ref.theta(:,:,f);
end
ref.theta = wrapToPi(ref.theta);

% bia, randomly distributed in [-.5,.5]
ref.bia = rand(npattern,nframe) - .5;

%% generate responds with constant shifting

% beta, follow cauchy distribution
ref.beta = ones(npattern,ntrans,nframe);

% theta, start from one frame and spread to all frames
% follow Gaussian distribution
ref.theta = repmat(reshape(linspace(0,2*pi,nframe),[1,1,nframe]),[npattern,ntrans,1]);

% bia, randomly distributed in [-.5,.5]
ref.bia = rand(npattern,nframe) - .5;

%% generate animation

npixel = prod(frmsz);
% reshape parameters for calculating
alpha = reshape(ref.alpha,[npixel,npattern,1,1]);
phi   = reshape(ref.phi,[npixel,1,ntrans,1]);
beta  = reshape(ref.beta,[1,npattern,ntrans,nframe]);
theta = reshape(ref.theta,[1,npattern,ntrans,nframe]);
bia   = reshape(ref.bia,[1,npattern,1,nframe]);

% initialize video parameters
v.res = frmsz;
v.ffindex = 1;
% generate video by generative model
v.v = reshape(genmodel(alpha,phi,beta,theta,bia),[npixel,nframe]);
% % add noise
% v.v = v.v + randn(npixel,nframe) * sigma.noise;

%% calculate objective value of reference model
delta = reshape(v.v,[npixel,1,1,nframe]) - genmodel(alpha,phi,beta,theta,bia);
ref.obj = objFunc(alpha,phi,beta,theta,bia,delta,sigma,v.ffindex,v.res);

clear alpha phi beta theta bia delta

%% ========================================================================
% train transform-pattern model to reconstruct reference model

%% initialize model from random variable
m = tpmodel(v,'nepoch',100,'nadapt',1,'ninfer',1,'npattern',npattern,'ntrans',ntrans, ...
    'noiseprior',sigma.noise,'sparseprior',sigma.sparse,'slowprior',sigma.slow, ...
    'patternbasesmoothprior',sigma.smpat,'transformbasesmoothprior',1, ...
    'NegativeCutProbability',0.9,'AlphaNormalization',false);

%% initialize model from reference with some noise
m = ref;
m.alpha = m.alpha + randn(size(m.alpha)) * mnoise;
m.phi   = m.phi + randn(size(m.phi)) * mnoise;
m.beta  = m.beta + randn(size(m.beta)) * mnoise;
m.theta = m.theta + randn(size(m.theta)) * mnoise;
m.bia   = m.bia + randn(size(m.bia)) * mnoise;

m.sigma = ref.sigma;

%% retrain model
[m,v] = tpmodel(v,'model',m,'nepoch',100,'nadapt',1,'ninfer',1, ...
    'NegativeCutProbability',.3,'AlphaNormalization',false);

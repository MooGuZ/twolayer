% This file aims at debuging the programs of transform-pattern model by
% creating sample bases and responds according to the distribution. This
% file is organized by cell units.
%
% MooGu Z. <hzhu@case.edu>
% Feb 19, 2015 - Version 0.1

clear all
close all

rng shuffle

%% fundamental settings
% number of bases
npattern = 1;
ntrans   = 1;
% size of frame
frmsz = [32,32];
% number of frams
nframe = 24;
% model noise
mnoise = 0.1;
% parameters of prior and likelihood
sigma.noise   = 0.1;
sigma.sparse  = 1;
sigma.slow    = 2*pi;
sigma.smpat   = 1;
sigma.smtrans = 2*pi;
% save sigma setting to reference model
ref.sigma = sigma;

%% generate reference bases
ref.alpha = reshape(pbaseGen([.5,.5],13,frmsz(1)),[prod(frmsz),1]);
ref.phi   = reshape(tbaseGen([1,1],3,frmsz(1)),[prod(frmsz),1]);

%% plot reference bases
imshow(baseplot(ref.alpha,ref.phi,frmsz));

%% generate responds randomly

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
ref.theta = reshape(linspace(0,2*pi,nframe),[1,1,nframe]);

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
video.res = frmsz;
video.ffindex = 1;
% generate video by generative model
video.v = reshape(genmodel(alpha,phi,beta,theta,bia),[npixel,nframe]);
% % add noise
% video.v = video.v + randn(npixel,nframe) * sigma.noise;

%% calculate objective value of reference model
delta = reshape(video.v,[npixel,1,1,nframe]) - genmodel(alpha,phi,beta,theta,bia);
ref.obj = objFunc(alpha,phi,beta,theta,bia,delta,sigma,video.ffindex,video.res);

%% ========================================================================
% train transform-pattern model to reconstruct reference model

%% initialize model from random variable
m = tpmodel(video,'nepoch',10,'nadapt',30,'ninfer',30,'npattern',npattern,'ntrans',ntrans, ...
    'noiseprior',sigma.noise,'sparseprior',sigma.sparse,'slowprior',sigma.slow, ...
    'patternbasesmoothprior',sigma.smpat,'transformbasesmoothprior',sigma.smpat);

%% initialize model from reference with some noise
m = ref;
m.alpha = m.alpha + randn(size(m.alpha)) * mnoise;
m.phi   = m.phi + randn(size(m.phi)) * mnoise;
m.beta  = m.beta + randn(size(m.beta)) * mnoise;
m.theta = m.theta + randn(size(m.theta)) * mnoise;
m.bia   = m.bia + randn(size(m.bia)) * mnoise;

%% retrain model
[m,v] = tpmodel(video,'model',m,'nepoch',100,'nadapt',70,'ninfer',70);

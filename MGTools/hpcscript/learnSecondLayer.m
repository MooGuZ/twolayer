% clear all; close all; clc;

% Add Pathes
root = [pwd,'/'];
addpath([root,'data'],[root,'state'],[root,'init']);
addpath([root,'code'], ...
    [root,'code/MGTools'],[root,'code/MGTools/io'],[root,'code/MGTools/learn'], ...
    [root,'code/MGTools/showrst'],[root,'code/MGTools/transform'], ...
    [root,'code/tools'],[root,'code/tools/minFunc_ind'], ...
    [root,'code/tools/minFunc_ind/logistic']);

% Load Initial State
load init/init.mat

% Close GPU option
p.use_gpu = true;

% Initial Random Number Generator
rng('shuffle');

warning('off','MATLAB:nearlySingularMatrix');

% Setup Parameters
param.epoches = [13];
param.tfactor = [1];
param.nsave = 10000;
param.testsz = 200;
param.rfile = './data/R1L-PATCH32-21-Apr-2014 12:06:50.mat';
% Learn Bases in First Layer
[m,p,snr] = learnPMCode(m,p,param);

% Reopen GPU option
% p.use_gpu = true;

% Rewrite Initial States
save init/init.mat m p
save state/snr-20140422.mat snr

% END
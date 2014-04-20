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

% Open GPU option
p.use_gpu = true;

% Learn Bases in First Layer
[m,p]=learnBases(m,p,16);

% Rewrite Initial States
save init/init.mat m p

% END
% This script file initialize Complex Bases Model and trying to learn complex bases
% fucntions from Fix View-Point dataset
%
% MooGu Z. <hzhu@case.edu>
% June 17, 2014 - Version 0.1

root = [pwd,'/'];
addpath([root,'data'],[root,'state'],[root,'init']);
addpath([root,'code'], ...
    [root,'code/MGTools'],[root,'code/MGTools/io'],[root,'code/MGTools/learn'], ...
    [root,'code/MGTools/init'],[root,'code/MGTools/showrst'],[root,'code/MGTools/transform'], ...
    [root,'code/tools'],[root,'code/tools/minFunc_ind'], ...
    [root,'code/tools/minFunc_ind/logistic']);

frmSize  = 32;
dataType = 'fvp';
dataPath = [root,'data/fvp-20140602'];

% Initialize Complex Bases Model
[m,p] = modelInitParameters(frmSize,dataType,dataPath);
if p.whitening.enable
    [Data,m,p] = modelInitWhitening(m,p);
end
[m,p] = modelInitBases(m,p);

% Utilize GPU's Power
p.use_gpu = true;

% Initialize Random Number Generator
rng('shuffle');

% Learn Firstlayer Bases Functions
[m,p] = learnCBases(m,p,48,10000);
p.firstlayer.eta_dA_target = .25 * p.firstlayer.eta_dA_target;
[m,p] = learnCBases(m,p,32,10000);

% Save Final Status
save init/final.mat m p

% END
function ana = anaMergeByModel(A,B)
% ANAMERGEBYDATA merge two analysis structure according to their identical
% model profiles.
%
% see also modelAnalysis, anaParamSetup, anaPlot, anaMergeByData, and
% anaUpdate.
%
% MooGu Z. <hzhu@case.edu>
% ------------------------------------------------------------------------
% Version 0.1 [June 30, 2014] - Start Version

% check completeness of analysis structures
if A.nproc ~= A.data.num
    disp(['The first analysis structure is not complete.\n', ...
        'Please complete it by MODELANALYSIS at first']);
elseif B.nproc ~= B.data.num
    disp(['The second analysis structure is not complete.\n', ...
        'Please complete it by MODELANALYSIS at first']);
end

% check consistency of dataset
if A.model.num ~= B.model.num || ~isempty(setxor(A.model.nameList,B.model.nameList))
    error('Model profiles are not consistent, merge operation aborted!');
end

% Basic information
ana.nproc    = A.nproc;
ana.fileName = 'anaMerged.mat';

% Data Information
ana.data = A.data;

% Check and set model path and namelist
if strcmp(A.model.path,B.model.path)
    ana.model.path = A.model.path;
    ana.model.nameList = [A.model.nameList,B.model.nameList];
else
    ana.model.nameList = [ ...
        cellfun(@(c)([A.model.path,c]),A.model.nameList,'UniformOutput',false), ...
        cellfun(@(c)([B.model.path,c]),B.model.nameList,'UniformOutput',flase)];
end

% Model Inforamtion
ana.model.IDList = [A.model.IDList,B.model.IDList];
ana.model.nbase = A.model.nbase;
ana.model.num = A.model.num + B.model.num;

% check whether or not model need reorganizing
swModelReorg = A.model.IDList(end) > B.model.IDList(1);

% if reorganization is necessary find the order
if swModelReorg
    [~,modelMap] = sort(ana.modelIDList,'ascend');
end

% Find Data Mapping Index
[~,sortIndex] = sort(B.data.index,'ascend');
dataMap = sortIndex(A.data.index);

% Noise Information
ana.noise.energy = [A.noise.energy; B.noise.energy(:,dataMap)];
ana.noise.likelihood = [A.noise.likelihood; B.noise.likelihood(:,dataMap)];
if swModelReorg
    ana.noise.energy = ana.noise.energy(modelMap,:);
    ana.noise.likelihood = ana.noise.likelihood(modelMap,:);
end

% Responds Distribution
ana.resp.dist.nbox = A.resp.dist.nbox;
ana.resp.dist.unit = A.resp.dist.unit;
ana.resp.dist.crd  = A.resp.dist.crd;
ana.resp.dist.bsz  = A.resp.dist.bsz;
ana.resp.dist.amp  = [A.resp.dist.amp; B.resp.dist.amp];
ana.resp.dist.damp = [A.resp.dist.damp; B.resp.dist.damp];
if swModelReorg
    ana.resp.dist.amp = ana.resp.dist.amp(modelMap,:);
    ana.resp.dist.damp = ana.resp.dist.damp(modelMap,:);
end

% Responds Probability Description
ana.resp.prob.sparse = [A.resp.prob.sparse; B.resp.prob.sparse(:,dataMap)];
ana.resp.prob.slow   = [A.resp.prob.slow; B.resp.prob.slow(:,dataMap)];
if swModelReorg
    ana.resp.prob.sparse = ana.resp.prob.sparse(modelMap,:);
    ana.resp.prob.slow = ana.resp.prob.slow(modelMap,:);
end

% Tophit List
if A.model.IDList(end) > B.model.IDList(end)
    ana.tophit = A.tophit;
else
    ana.tophit = B.tophit;
end

% Generate Summary
ana.summary.SNR = -10*log10(sum(ana.noise.energy(:,1:ana.nproc),2) ...
    / sum(ana.data.energy(1:ana.nproc)));
ana.summary.likelihood = mean(ana.noise.likelihood(:,1:ana.nproc),2);
ana.summary.sparse = mean(ana.resp.prob.sparse(:,1:ana.nproc),2);
ana.summary.slow = mean(ana.resp.prob.slow(:,1:ana.nproc),2);
ana.summary.prob = ana.summary.sparse + ana.summary.slow ...
    + ana.summary.likelihood;
% Save Result
ana.finishtime = datestr(now);
save(ana.fileName,'ana');

end
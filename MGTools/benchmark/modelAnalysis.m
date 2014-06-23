function ana = modelAnalysis(inArgA,inArgB)
% MODELANALYSIS produce evaluations on the basis functions learned in complex bases
% model. The analysis process would based on the dataset provided in parameters.
%
%   ana = modelAnalysis(modelp,datap) initialize analysis process by model
%   parameters MODELP and data parameters DATAP. ANA is the structure contains all
%   the information about the evaluation.
%
%   ana = modelAnalysis(ana) continue and finish an analysis process.
%
% see also anaParamSetup, anaPlot.
%
% MooGu Z. <hzhu@case.edu>
% ----------------------------------------------------------------------------------
% Version 0.1 [June 18, 2014] - Start Version
% Version 0.2 [June 21, 2014] - Change input argumen into general form 
% Version 0.3 [June 23, 2014] - Add support to log-scale distribution plots

% if no ANA exist, initialize one
switch (nargin)
  case 1
    ana = inArgA;
    fprintf('Work Start @ %s with %d Video Clips have been Analyzed.\n', ...
        datestr(now),ana.nproc);
  case 2
    ana.nproc = 0;
    ana.fileName = [inArgB.anaFilePath,'analysis-',datestr(now),'.mat'];
    
    ana.data.path = inArgB.path;
    ana.data.nameList = inArgB.nameList;
    ana.data.batchSize = inArgB.batchSize;
    ana.data.npixel = inArgB.npixel;
    ana.data.nframe = inArgB.nframe;
    
    ana.data.num = numel(inArgB.nameList);
    ana.data.index = randperm(ana.data.num);
    ana.data.energy = zeros(1,ana.data.num);
    
    ana.model.path = inArgA.path;
    ana.model.nameList = inArgA.nameList;
    ana.model.IDList = inArgA.IDList;
    ana.model.nbase = inArgA.baseQuantity;
    
    ana.model.num = numel(inArgA.nameList);
    
    ana.noise.energy = zeros(ana.model.num,ana.data.num);
    ana.noise.likelihood = zeros(ana.model.num,ana.data.num);
    
    ana.resp.dist.nbox = inArgA.distBoxNum;
    ana.resp.dist.unit = inArgA.unitOrder;

    % calculate coordinates and box size of distribution plots
    ana.resp.dist.crd = logspace(ana.resp.dist.unit,1,ana.resp.dist.nbox);    
    ana.resp.dist.bsz = [ ...
        .5 * sum(ana.resp.dist.crd(1:2)), ...
        diff(conv(ana.resp.dist.crd,[.5,.5],'valid'),1,2), ...
        Inf];
    
    ana.resp.dist.amp = zeros(ana.model.num,ana.resp.dist.nbox);
    ana.resp.dist.damp = zeros(ana.model.num,ana.resp.dist.nbox);
    ana.resp.prob.sparse = zeros(ana.model.num,ana.data.num);
    ana.resp.prob.slow = zeros(ana.model.num,ana.data.num);
    
    ana.tophit = zeros(ana.model.nbase,ana.data.num);
    
    fprintf('Work Start @ %s\n',datestr(now));
  otherwise
    error('Wrong Syntax in using function MODELANALYSIS.');  
end

% Create Shortcuts
ndata   = ana.data.num;
nmodel  = ana.model.num;
dpath   = ana.data.path;
dindex  = ana.data.index;

while(ana.nproc < ndata)
    nbatch = min(ana.data.batchSize, ndata-ana.nproc);
    % load the data that need to be processed in this batch
    Data = zeros(ana.data.npixel,ana.data.nframe,nbatch);
    for j = 1 : nbatch
        Data(:,:,j) = gif2anim( ...
            [dpath,ana.data.nameList{dindex(ana.nproc+j)}], ...
            1:ana.data.nframe);
    end
    % apply each model profile on these data
    for i = 1 : nmodel
        % load model profile
        clear m p
        load([ana.model.path,ana.model.nameList{i}]);
        % Shortcuts
        beta  = p.firstlayer.a_cauchy_beta;
        sigma = p.firstlayer.a_cauchy_sigma;
        lambda = p.firstlayer.a_lambda_S;
        % Whiten the Data
        if i == 1
            Data = reshape(Data,ana.data.npixel,ana.data.nframe*nbatch);
            Data = m.whitenMatrix * bsxfun(@minus,Data,m.imageMean);
            Data = reshape(Data,m.M,ana.data.nframe,nbatch);
            for j = 1 : nbatch
                ana.data.energy(ana.nproc+j) = sum(sum(Data(:,:,j).^2));
            end
        end
        % ======= modify model profile =======
        p.use_gpu = exist('gsingle','builtin');
        % ====================================
        for j = 1 : nbatch
            [Z,noise] = infer_Z(Data(:,:,j),m,p);
            
            ana.noise.energy(i,ana.nproc+j) = sum(noise(:).^2);
            ana.noise.likelihood(i,ana.nproc+j) = ...
                sum(sum(bsxfun(@times,.5*m.I_noise_factors,noise.^2))); 
            
            amp = abs(Z); damp = diff(amp,1,2);
            
            ana.resp.dist.amp(i,:) = ana.resp.dist.amp(i,:) ...
                + hist(amp(:),ana.resp.dist.crd);
            ana.resp.dist.damp(i,:) = ana.resp.dist.damp(i,:) ...
                + hist(abs(damp(:)),ana.resp.dist.crd);
            
            ana.resp.prob.sparse(i,ana.nproc+j) = ...
                beta * sum(log(1+(amp(:)/sigma).^2));
            ana.resp.prob.slow(i,ana.nproc+j) = ...
                .5 * lambda * sum(damp(:).^2);
            
            if i == nmodel
                [~,ana.tophit(:,ana.nproc+j)] = sort(sum(amp,2),'descend');
            end            
        end
    end
    ana.nproc = ana.nproc + nbatch;
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
    fprintf('%4d Video Clips Have Been Analyzed @ %s\n', ...
        ana.nproc,datestr(now));
end
    
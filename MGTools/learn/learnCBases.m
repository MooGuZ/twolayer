function [m,p] = learnCBases(m,p,nEpoch,nSave)
% LEARNBASES learn bases in first layer (This is a version for Fix View-Point
% Dataset)
%
% [m,p] = learnCBases(m,p,nEpoch,nSave)
%
% MooGu Z. <hzhu@case.edu>
% June 17, 2014 - Version 0.1

% number of iteration per save
if ~exist('nSave','var')
    nSave = 10000;
end

% Prepare for GPU computing
if p.use_gpu
    m.A = gsingle(m.A);
end

% Load Data into Memory Only Once
if p.data.scope >= p.data.quantity
    fprintf('Loading %d video clips ... ', p.data.scope);
    Data = loadDataBatch(m,p);
    fprintf('Done @ %s\n', datestr(now()));
end

fprintf('\nFirst-Layer-Bases Learning (%2d Epoches) start @ %s\n',nEpoch,datestr(now));
for iEpoch = 1 : nEpoch
    % Load Data for each Epoch, if necessary
    if p.data.scope < p.data.quantity
        fprintf('Loading %d video clips ... ', p.data.scope);
        Data = loadDataBatch(m,p);
        fprintf('Done @ %s\n', datestr(now()));
    end

    fprintf('Learning '); infotag = 0.1;
    % Generate random sequence
    index = randpermext(p.data.scope,nSave);
    for iter = 1 : nSave
        if p.use_gpu
            video = gsingle(im2double(Data(:,(index(iter)-1)*p.data.nframe+1 : ...
                                 index(iter)*p.data.nframe)));
        else
            video = im2double(Data(:,(index(iter)-1)*p.data.nframe+1 : ...
                                   index(iter)*p.data.nframe));
        end
        % Apply Whitening Operation
        if p.whitening.enable
            video = m.whitenMatrix * bsxfun(@minus,video,m.imageMean);
        end
        % Infer Complex Bases Responds
        [Z,I_E] = infer_Z(video,m,p);
        % Adapt Complex Bases
        [m,p] = adapt_firstlayer(Z,I_E,m,p);
        m.t(1) = m.t(1) + 1;
        % update console output
        if iter/nSave >= infotag
            infotag = infotag + 0.1;
            fprintf('.');
        end
    end
    iterstr = strrep(mat2str(m.t),' ',',');
    % save learning states
    save_model([p.autosave.path,'FirstLayerBases-Iteration', ...
                iterstr,'.mat'],m,p);
    % show information
    fprintf(' Iteration %s DONE @ %s\n', iterstr, datestr(now()));
end

end
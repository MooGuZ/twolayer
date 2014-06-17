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
    Data = loadDataBatch(p.data);
end

fprintf('\nFirst-Layer-Bases Learning (%2d Epoches) start @ %s\n',nEpoch,datestr(now));
for iEpoch = 1 : nEpoch
    % Load Data for each Epoch, if necessary
    if p.data.scope < p.data.quantity
        Data = loadDataBatch(p.data);
    end
    % Generate random sequence
    index = randi(p.data.scope,nSave,1);
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
        [Z,I_E,fexit] = infer_Z(video,m,p);
        % Adapt Complex Bases
        [m,p] = adapt_firstlayer(Z,I_E,m,p);
        m.t(1) = m.t(1) + 1;
    end
    iterstr = strrep(mat2str(m.t),' ',',');
    % save learning states
    save_model([p.autosave.path,'state/FirstLayerBases-Iteration',iterstr,'.mat'],m,p);
    % Output Infomation in Console
    disp(['Learning Iteration ',iterstr,' DONE @ ',datestr(now)]);
end

end
function [m,p] = learnBases(m,p,nEpoch,nSave)
% LEARNBASES learn bases in first layer

% number of iteration per save
if ~exist('nSave','var')
    nSave = 10000;
end
% number of chunks read in one time
nChunks = ceil(nSave / p.patches_load);

% Prepare for GPU computing
if p.use_gpu
    m.A = gsingle(m.A);
end

% Parameters for Data Frames
vsize  = p.imsz-2*p.BUFF-p.topmargin;
hsize  = p.imsz-2*p.BUFF;
vpos   = 1 + p.topmargin+p.BUFF;
hpos   = 1 + p.BUFF;
vrange = vsize - m.patch_sz + 1;
hrange = hsize - m.patch_sz + 1;

% Read chunks into memory in advance
chunks = cell(p.num_chunks,1);
for i = 1 : p.num_chunks
    chunks{i} = readdata(m,p,'chunks',i,'patchsize',[vsize,hsize],'position',[vpos,hpos]);
end

fprintf('\nFirst-Layer-Bases Learning (%2d Epoches) start @ %s\n',nEpoch,datestr(now));
for iEpoch = 1 : nEpoch
    % Generate random chunk sequence
    index = randi(p.num_chunks,nChunks,1);
    for c = 1 : nChunks
        % Put chunk into GPU's memory
        if p.use_gpu
            D = gsingle(reshape(chunks{index(c)},vsize,hsize,p.imszt));
        else
            D = reshape(chunks{index(c)},vsize,hsize,p.imszt);
        end
        % Learn Bases from current chunk
        for iter = 1 : p.patches_load
            fexit = false;
            while ~fexit
                vs = ceil(vrange*rand(1));
                hs = ceil(hrange*rand(1));
                X  = reshape(D(vs:vs+m.patch_sz-1,hs:hs+m.patch_sz-1,:), ...
                    m.patch_sz^2,p.imszt);
                X = m.whitenMatrix * bsxfun(@minus,X,m.imageMean);
                [Z,I_E,fexit] = infer_Z(X,m,p);
            end
            [m,p] = adapt_firstlayer(Z,I_E,m,p);
            m.t(1) = m.t(1) + 1;
        end
    end
    iterstr = strrep(mat2str(m.t),' ',',');
    % save learning states
    save_model([p.autosave_path,'FirstLayerBases-Iteration',iterstr,'.mat'],m,p);
    % Output Infomation in Console
    disp(['Learning Iteration ',iterstr,' DONE @ ',datestr(now)]);
end

end
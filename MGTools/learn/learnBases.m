function [m,p] = learnBases(m,p,nEpoch)
% LEARNBASES learn bases in first layer

% number of iteration per save
nSave = 5000;
% number of chunks read in one time
nChunks = ceil(nSave / p.patches_load);

% Prepare for GPU computing
if p.use_gpu
    m.A = gsingle(m.A);
end

% Parameters for Data Frames
vsize  = p.imsz-2*p.BUFF-p.topmargin;
hsize  = p.imsz-2*p.BUFF;
vpos   = 1 + p.topmarin+p.BUFF;
hpos   = 1 + p.BUFF;
vrange = vsize - m.patch_sz + 1;
hrange = hsize - m.pathc_sz + 1;

for iEpoch = 1 : nEpoch
    % Read Data
    data  = readdata(m,p,'chunks',ceil(rand(nChunks,1)*p.num_chunks),...
                     'patchsize',[vsize,hsize],'position',[vpos,hpos]);
    for c = 1 : nChunks
        % Put chunk into GPU's memory
        if p.use_gpu
            D = gsingle(data(:,:,(c-1)*p.imszt+1:c*p.imszt));
        end
        % Learn Bases from current chunk
        for iter = 1 : p.patch_load
            fexit = false
            while ~fexit
                vs = ceil(vrange*rand(1));
                hs = ceil(hrange*rand(1));
                [Z,I_E,fexit] = ...
                    infer_Z(D(vs:vs+m.patch_sz-1,hs:hs+m.patch_sz-1,:),m,p);
            end
            [m,p] = adapt_firstlayer(Z,I_E,m,p);
            m.t = m.t + 1;
        end
    end
    % save learning states
    save_model([p.autosave_path,'Iteration(',num2str(m.t),').mat'],m,p);
end
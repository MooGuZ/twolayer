function [m,p] = learnBases(m,p,nEpoch)
% LEARNBASES learn bases in first layer

% number of iteration per save
nSave = 20;
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

disp(['First-Layer-Bases Learning (',num2str(nEpoch),' Epoches) start @ ',datestr(now)]);
for iEpoch = 1 : nEpoch
    % Read Data
    data  = readdata(m,p,'chunks',ceil(rand(nChunks,1)*p.num_chunks),...
        'patchsize',[vsize,hsize],'position',[vpos,hpos]);
    for c = 1 : nChunks
        % Put chunk into GPU's memory
        if p.use_gpu
            D = gsingle(reshape(data(:,(c-1)*p.imszt+1:c*p.imszt),vsize,hsize));
        else
            D = reshape(data(:,(c-1)*p.imszt+1:c*p.imszt),vsize,hsize,p.imszt);
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
            m.t = m.t + 1;
        end
    end
    % save learning states
    save_model([p.autosave_path,'Iteration(',num2str(m.t),').mat'],m,p);
    % Output Infomation in Console
    disp(['Iteration(',num2str(m.t),') DONE @ ',datestr(now)]);
end
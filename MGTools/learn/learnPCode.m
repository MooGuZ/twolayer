function [m,p] = learnPCode(m,p,nEpoch,nSave)
% LEARNPCODE learn pattern and motion bases (codes) in second layer

% number of iteration per save
if ~exist('nSave','var'), nSave = 10000; end

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

fprintf('\nPattern-Bases Learning (%2d Epoches) start @ %s\n',nEpoch,datestr(now));
fprintf('Generating First Layer Responds ... ')

% number of chunks read into memory
nSegments = ceil(p.load_segments / p.num_chunks);
% Initialize storage for Amplitude and Phase
logAmp = zeros(m.N,p.imszt*nSegments*p.num_chunks);
dPhase = Amp;
for c = 1 : p.num_chunks
    % Put chunk into GPU's memory
    if p.use_gpu
        D = gsingle(reshape(chunks{index(c)},vsize,hsize,p.imszt));
    else
        D = reshape(chunks{index(c)},vsize,hsize,p.imszt);
    end
    % Calculate Responds from First Layer
    for i = 1 : nSegments
        fexit = false;
        while ~fexit
            vs = ceil(vrange*rand(1));
            hs = ceil(hrange*rand(1));
            X  = reshape(D(vs:vs+m.patch_sz-1,hs:hs+m.patch_sz-1,:), ...
                m.patch_sz^2,p.imszt);
            X = m.whitenMatrix * bsxfun(@minus,X,m.imageMean);
            [Z,~,fexit] = infer_Z(X,m,p);
        end
        sind = ((c-1)*nSegments+i-1)*p.imszt + 1;
        eind = ((c-1)*nSegments+i)*p.imszt;
        logAmp(:,sind:eind) = single(log(abs(Z)));
        dPhase(:,sind:eind) = single([angle(Z(:,1)),angle(Z(:,2:end))-angle(Z(:,1:end-1))]);
    end
end

% Release Chunks
clear('chunks');

% Get statistic profiles of amplitude
m.loga_means = zeros(m.N,1);
m.loga_factors = m.loga_means;
for i = 1 : m.N
    m.loga_means(i) = mean(logAmp(i,:));
    m.loga_factors(i) = sqrt(1/(10*var(logAmp(i,:))));
end

end
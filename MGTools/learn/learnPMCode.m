function [m,p] = learnPMCode(m,p,nEpoch,nSave,respfile)
% LEARNPCODE learn pattern and motion bases (codes) in second layer

% number of iteration per save
if ~exist('nSave','var'), nSave = 10000; end

% Prepare for GPU computing
if p.use_gpu
    m.A = gsingle(m.A);
end

if exist('respfile','var')
    load(respfile);
else
    fprintf('Generating First Layer Responds start @ %s\n',datestr(now));
    % Parameters for Data Frames
    vsize  = p.imsz-2*p.BUFF-p.topmargin;
    hsize  = p.imsz-2*p.BUFF;
    vpos   = 1 + p.topmargin+p.BUFF;
    hpos   = 1 + p.BUFF;
    vrange = vsize - m.patch_sz + 1;
    hrange = hsize - m.patch_sz + 1;
    
    % Read chunks into memory in advance
    fprintf('Reading Data Chunks ... ');
    chunks = cell(p.num_chunks,1);
    for i = 1 : p.num_chunks
        chunks{i} = readdata(m,p,'chunks',i,'patchsize', ...
            [vsize,hsize],'position',[vpos,hpos]);
    end
    fprintf('DONE\n');
    
    % number of chunks read into memory
    nSegments = ceil(p.load_segments / p.num_chunks);
    p.load_segments = nSegments * p.num_chunks;
    % Initialize storage for Amplitude and Phase
    R1L.logAmp = zeros(m.N,p.imszt*nSegments*p.num_chunks,'single');
    R1L.dPhase = zeros(m.N,p.imszt*nSegments*p.num_chunks,'single');
    for c = 1 : p.num_chunks
        fprintf('Chunk %2d ... ',c);
        % Put chunk into GPU's memory
        if p.use_gpu
            D = gsingle(reshape(chunks{c},vsize,hsize,p.imszt));
        else
            D = single(reshape(chunks{c},vsize,hsize,p.imszt));
        end
        % Calculate Responds from First Layer
        for i = 1 : nSegments
            vs = ceil(vrange*rand(1));
            hs = ceil(hrange*rand(1));
            X  = reshape(D(vs:vs+m.patch_sz-1,hs:hs+m.patch_sz-1,:), ...
                m.patch_sz^2,p.imszt);
            X  = m.whitenMatrix * bsxfun(@minus,X,m.imageMean);
            Z  = infer_Z(X,m,p);
            sind = ((c-1)*nSegments+i-1)*p.imszt + 1;
            eind = ((c-1)*nSegments+i)*p.imszt;
            R1L.logAmp(:,sind:eind) = log(abs(Z));
            R1L.dPhase(:,sind:eind) = ...
                [angle(Z(:,1)),angle(Z(:,2:end))-angle(Z(:,1:end-1))];
        end
        fprintf('DONE (%d Segments Loaded)\n',i);
    end
    
    % Release Chunks
    clear('chunks');
    
    % Get statistic profiles of amplitude
    m.loga_means = zeros(m.N,1);
    m.loga_factors = m.loga_means;
    % Calculating in For-Loop for memory issues
    for i = 1 : m.N
        m.loga_means(i) = mean(R1L.logAmp(i,:));
        m.loga_factors(i) = sqrt(1/(10*var(R1L.logAmp(i,:))));
    end
    
    % Save First Layer Responds
    respfile = ['R1L-PATCH',num2str(m.patch_sz),'-',datestr(now),'.mat'];
    save(['data/',respfile],'R1L','m','p','-v7.3');
    fprintf('Save responds of first layer into %s\n',respfile);
end

fprintf('\nPattern-Bases Learning (%2d Epoches) start @ %s\n',nEpoch,datestr(now));


fprintf('Pattern-Bases Learning Process DONE @ %s\n',datestr(now));
end

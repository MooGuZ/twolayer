function [m,p,respfile] = cbaserespond(m,p,nSegments)
% CBASERESPOND collect complex bases responds for movies
%
%   [m,p,respfile] = cbaserespond(m,p,nSegments)
%
% MooGu Z. <hzhu@case.edu>
% April 20, 2014 - Version 0.1


% Number of segments loaded for each data chunk
if ~exist('nSegments','var'), nSegments = 130; end
p.load_segments = nSegments * p.num_chunks;

% Prepare for GPU computing
if p.use_gpu
    m.A = gsingle(m.A);
end

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

% Initialize storage for Amplitude and Phase
R1L.logAmp = zeros(m.N,p.imszt*nSegments*p.num_chunks,'single');
R1L.dPhase = zeros(m.N,p.imszt*nSegments*p.num_chunks,'single');
for c = 1 : p.num_chunks
    fprintf('Chunk %2d ... ',c);
    % Put chunk into memory
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
        R1L.logAmp(:,sind:eind) = log(abs(Z)+eps);
        R1L.dPhase(:,sind:eind) = ...
            [angle(Z(:,1)),angle(Z(:,2:end))-angle(Z(:,1:end-1))];
    end
    fprintf('DONE (%d Segments Loaded)\n',i);
end

% Release Chunks
clear('chunks');

% Get statistic profiles of amplitude
m.loga_means = mean(R1L.logAmp,2);
m.loga_factors = sqrt(0.1./var(R1L.logAmp,0,2));

% Type Transformation
if p.use_gpu
    m.A = double(m.A);
end
% Save First Layer Responds
respfile = ['R1L-PATCH',num2str(m.patch_sz),'-',datestr(now),'.mat'];
save(['data/',respfile],'R1L','m','p','-v7.3');
fprintf('Save responds of first layer into %s @ %s\n',respfile,datestr(now));
respfile = ['./data/',respfile];

end

function [m,p] = collectFirstLayerRespond(m,p,respfile)
% CBASERESPOND collect complex bases responds for movies
%
% MooGu Z. <hzhu@case.edu>
% Oct 12, 2015


% % Number of segments loaded for each data chunk
% if ~exist('nSegments','var'), nSegments = 130; end
% p.load_segments = nSegments * p.num_chunks;

% Prepare for GPU computing
if p.use_gpu
    m.A = gsingle(m.A);
end

fprintf('Generating First Layer Responds start @ %s\n',datestr(now));
% Initialize storage for Amplitude and Phase
R1L.logAmp = zeros(m.N,p.data.nframe*p.data.quantity,'single');
R1L.dPhase = zeros(m.N,p.data.nframe*p.data.quantity,'single');
% inference one by one
fprintf('INFER PROCESS '); infotag = 0.1;
for i = 1 : p.data.quantity
    X = im2single(gif2anim(fullfile(p.data.path, p.data.nameList{i}), 1 : p.data.nframe));
    X = m.whitenMatrix * bsxfun(@minus, X, m.imageMean);
    Z = infer_Z(X, m, p);
    ind = (i - 1) * p.data.nframe + 1 : i * p.data.nframe;
    R1L.logAmp(:, ind) = log(abs(Z) + eps);
    R1L.dPhase(:, ind) = [angle(Z(:, 1)), ...
        wrapToPi(angle(Z(:, 2:end)) - angle(Z(:, 1:(end-1))))];
    if (i / p.data.quantity) >= infotag
        infotag = infotag + 0.1;
        fprintf('>');
    end
end
fprintf(' DONE');

% Get statistic profiles of amplitude
m.loga_means = mean(R1L.logAmp,2);
m.loga_factors = sqrt(0.1./var(R1L.logAmp,0,2));

% Type Transformation
if p.use_gpu
    m.A = double(m.A);
end
% Save First Layer Responds
save(respfile, 'R1L', 'm', 'p', '-v7.3');
fprintf('Save responds of first layer into %s @ %s\n',respfile,datestr(now));

end

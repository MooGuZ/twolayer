function [m,p,loglh,snr] = learnPMBase(m,p,nEpoch,param)
% LEARNPMCODE learn pattern and motion bases in second layer
%
%   [m,p,loglh,snr] = learnPMBase(m,p,nEpoch,param)
%
%   param
%    |- nsave   : number of iterations that save model paramter once
%    |- rfile   : file name of first layer responds of movie data
%    |- testsz  : size of test set
%
% MooGu Z. <hzhu@case.edu>
% May 28, 2014 - Version 0.1

% number of iteration per save
if isfield(param,'nsave') && ~isempty(param.nsave)
    nSave = param.nsave;
else
    nSave = 1000;
end

% Load First Layer Responds
if isfield(param,'rfile')
    m.rfile = param.rfile;
else
    [m,p,m.rfile] = cbaserespond(m,p);
end
load(m.rfile,'R1L');

% Initialized Performance Records
snr.logamp = zeros(nEpoch,1);
snr.dphase = zeros(nEpoch,1);
loglh.logamp = zeros(nEpoch,1);
loglh.dphase = zeros(nEpoch,1);

% Choose Test Set
if isfield(param,'testsz')
    testsz = param.testsz;
else
    testsz = 130;
end
testset = randi(p.load_segments,testsz,1);

% Segment Index Function
segs = @(n) (n-1)*p.imszt+1;
sege = @(n) n*p.imszt;

% Create Map of 1st-layer Bases coordinates in space and spacial frequency
if ~isfield(m,'Acoords'), m = fit_Acoords(m); end

fprintf('\nPattern-Bases Learning (%d Epoches) START @ %s\n',nEpoch,datestr(now));

% Convert Bases into GSINGLE
if p.use_gpu
    m.B = gsingle(m.B);
    m.D = gsingle(m.D);
end

for i = 1 : nEpoch
    I = randi(p.load_segments,nSave,1);
    for j = 1 : nSave
        if p.use_gpu
            P = gsingle(R1L.dPhase(:,segs(I(j))+1:sege(I(j))));
            A = gsingle(R1L.logAmp(:,segs(I(j)):sege(I(j))));
        else
            P = single(R1L.dPhase(:,segs(I(j))+1:sege(I(j))));
            A = single(R1L.logAmp(:,segs(I(j)):sege(I(j))));
        end
        % ==INFER MOTION CODES==
        % Roll Phase Difference into [-pi,pi]
        P(P < -pi) = P(P < -pi) + 2*pi;
        P(P >  pi) = P(P >  pi) - 2*pi;
        % Calculate Amplitude Mask
        mask = (A >= p.phasetrans.a_thresh);
        mask = mask(1:end-1) & mask(2:end);
        % Infering and Adapting
        W = inferMCode(P,mask,m,p);
        m = adaptMBase(W,P,mask,m,p);
        % Update counter
        m.t(2) = m.t(2) + 1;
        % ==INFER PATTERN CODES==
        % Centralize and Scale
        A = bsxfun(@minus,A,m.loga_means);
        A = bsxfun(@times,A,m.loga_factors);
        % Infering and Adapting
        V = inferPCode(A,m,p);
        m = adaptPBase(V,A,m,p);
        % Update counter
        m.t(3) = m.t(3) + 1;
    end
    % SNR Calculation
    for j = 1 : testsz
        % ==COLLECT DATA==
        if p.use_gpu
            P = gsingle(R1L.dPhase(:, ...
                segs(testset(j))+1:sege(testset(j))));
            A = gsingle(R1L.logAmp(:, ...
                segs(testset(j)):sege(testset(j))));
        else
            P = single(R1L.dPhase(:, ...
                segs(testset(j))+1:sege(testset(j))));
            A = single(R1L.logAmp(:, ...
                segs(testset(j)):sege(testset(j))));
        end
        % ==LEARN MOTION CODE BASES==
        % Roll Phase Difference into [-pi,pi]
        P(P < -pi) = P(P < -pi) + 2*pi;
        P(P >  pi) = P(P >  pi) - 2*pi;
        % Calculate Amplitude Mask
        mask = (A >= p.phasetrans.a_thresh);
        mask = mask(1:end-1) & mask(2:end);
        % Infering and Calculating SNR
        [W,likelihood] = inferMCode(P,mask,m,p);
        noise = sum((mask.*(P - m.D*W)).^2);
        snr.dphase(i) = snr.dphase(i) + ...
            10*log10(sum(P(:).^2)/noise);
        loglh.dphase(i) = loglh.dphase(i) + likelihood;
        % ==LEARN FORM CODE BASES==
        % Centralize and Scale
        A = bsxfun(@minus,A,m.loga_means);
        A = bsxfun(@times,A,m.loga_factors);
        % Infering and Calculating SNR
        [V,likelihood] = inferPCode(A,m,p);
        noise = evalPBase(m.B,V,P);
        snr.logamp(i) = snr.logamp(i) + ...
            10*log10(sum(A(:).^2)/noise);
        loglh.logamp(i) = loglh.logamp(i) + likelihood;
    end
    snr.dphase(i) = snr.dphase(i) / testsz;
    snr.logamp(i) = snr.logamp(i) / testsz;
    loglh.dphase(i) = loglh.dphase(i) / testsz;
    loglh.logamp(i) = loglh.logamp(i) / testsz;
    % Save Model and Parameters
    iterstr = strrep(mat2str(m.t),' ',',');
    % save learning states
    save_model([p.autosave_path,'state/PMCodeBases-Iteration', ...
        iterstr,'.mat'],m,p);
    % Output Infomation in Console
    fprintf('Learning Iteration %s DONE @ %s | SNR(A:%.2e,P:%.2e)\n', ...
        iterstr,datestr(now),snr.logamp(iepoch),snr.dphase(iepoch));
end

% Convert Bases into GSINGLE
if p.use_gpu
    m.B = double(m.B);
    m.D = double(m.D);
end

fprintf('Pattern-Bases Learning (%d Epoches) END @ %s\n', ...
    sum(nEpoch),datestr(now));

end

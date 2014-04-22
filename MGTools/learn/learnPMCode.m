function [m,p,snr] = learnPMCode(m,p,param)
% LEARNPMCODE learn pattern and motion bases (codes) in second layer
%
%   [m,p,snr] = learnPMCode(m,p,param)
%
%   param
%    |- epoches : array of numbers represent # of epochs for each procedure
%    |- tfactor : array of numbers corresponds scale factor of target error
%    |- nsave   : number of iterations that save model paramter once
%    |- rfile   : file name of first layer responds of movie data
%    |- testsz  : size of test set
%
% MooGu Z. <hzhu@case.edu>
% April 20, 2014 - Version 0.1

% number of iteration per save
if isfield(param,'nsave')
    nSave = param.nsave; 
else
    nSave = 10000;
end

% Load First Layer Responds
if isfield(param,'rfile')
    m.rfile = param.rfile;
else
    [m,p,m.rfile] = cbaserespond(m,p);
end
load(m.rfile,'R1L');

assert(all(size(param.epoches)==size(param.tfactor)), ...
    'ERROR : param.epoches need has same size of param.tfactor');
nEpoch  = param.epoches;
tFactor = param.tfactor;

% Initialize SNR
snr.logamp = zeros(sum(nEpoch),1);
snr.dphase = zeros(sum(nEpoch),1);

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

% Create Map of 1st-layer Bases coordinates in space and frequency
if ~isfield(m,'Acoords'), m = fit_Acoords(m); end

fprintf('\nPattern-Bases Learning (%d Epoches) START @ %s\n', ...
    sum(nEpoch),datestr(now));

% Convert Bases into GSINGLE
if p.use_gpu
    m.B = gsingle(m.B);
    m.D = gsingle(m.D);
end

for e = 1 : numel(nEpoch)
    % Scale Error Target in Learning
    p.phasetrans.eta_dD_target = p.phasetrans.eta_dD_target * tFactor(e);
    p.ampmodel.eta_dB_target = p.ampmodel.eta_dB_target * tFactor(e);
    for i = 1 : nEpoch(e)
        I = randi(p.load_segments,nSave,1);
        for j = 1 : nSave
            % ==COLLECT DATA==
            if p.use_gpu
                P = gsingle(R1L.dPhase(:,segs(I(j))+1:sege(I(j))));
                A = gsingle(R1L.logAmp(:,segs(I(j)):sege(I(j))));
                valid = gsingle(single(1));
            else
                P = single(R1L.dPhase(:,segs(I(j))+1:sege(I(j))));
                A = single(R1L.logAmp(:,segs(I(j)):sege(I(j))));
                valid = single(1);
            end
            % ==LEARN MOTION CODE BASES==
            % Roll Phase Difference into [-pi,pi]
            P(P < -pi) = P(P < -pi) + 2*pi;
            P(P >  pi) = P(P >  pi) - 2*pi;
            % Infering and Adapting
            [W,error] = infer_w(P,valid,m,p);
            [m,p] = adapt_phasetrans(W,error,m,p);
            % Update counter
            m.t(2) = m.t(2) + 1;
            % ==LEARN FORM CODE BASES==
            % Centralize and Scale
            A = bsxfun(@minus,A,m.loga_means);
            A = bsxfun(@times,A,m.loga_factors);
            % Infering and Adapting
            [V,error] = infer_v(A,m,p);
            [m,p] = adapt_ampmodel(V,[],error,m,p); 
            % Update counter
            m.t(3) = m.t(3) + 1;
        end
        % SNR Calculation
        iepoch = sum(nEpoch(1:e-1)) + i;
        for j = 1 : testsz
            % ==COLLECT DATA==
            if p.use_gpu
                P = gsingle(R1L.dPhase(:, ...
                    segs(testset(j))+1:sege(testset(j))));
                A = gsingle(R1L.logAmp(:, ...
                    segs(testset(j)):sege(testset(j))));
                valid = gsingle(1);
            else
                P = single(R1L.dPhase(:, ...
                    segs(testset(j))+1:sege(testset(j))));
                A = single(R1L.logAmp(:, ...
                    segs(testset(j)):sege(testset(j))));
                valid = single(1);
            end
            % ==LEARN MOTION CODE BASES==
            % Roll Phase Difference into [-pi,pi]
            P(P < -pi) = P(P < -pi) + 2*pi;
            P(P >  pi) = P(P >  pi) - 2*pi;
            % Infering and Calculating SNR
            [~,error] = infer_w(P,valid,m,p);
            snr.dphase(iepoch) = snr.dphase(iepoch) + ...
                10*log10(sum(P(:).^2)/sum(error(:).^2));
            % ==LEARN FORM CODE BASES==
            % Centralize and Scale
            A = bsxfun(@minus,A,m.loga_means);
            A = bsxfun(@times,A,m.loga_factors);
            % Infering and Calculating SNR
            [~,error] = infer_v(A,m,p);
            snr.logamp(iepoch) = snr.logamp(iepoch) + ...
                10*log10(sum(A(:).^2)/sum(error(:).^2));
        end
        snr.dphase(iepoch) = snr.dphase(iepoch) / testsz;
        snr.logamp(iepoch) = snr.logamp(iepoch) / testsz;
        % Save Model and Parameters
        iterstr = strrep(mat2str(m.t),' ',',');
        % save learning states
        save_model([p.autosave_path,'state/PMCodeBases-Iteration', ...
            iterstr,'.mat'],m,p);
        % Output Infomation in Console
        fprintf('Learning Iteration %s DONE @ %s | SNR(A:%.2e,P:%.2e)\n', ...
            iterstr,datestr(now),snr.logamp(iepoch),snr.dphase(iepoch));
    end
end

% Convert Bases into GSINGLE
if p.use_gpu
    m.B = double(m.B);
    m.D = double(m.D);
end

fprintf('Pattern-Bases Learning (%d Epoches) END @ %s\n', ...
    sum(nEpoch),datestr(now));

end

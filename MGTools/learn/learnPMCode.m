function [m,p,snr] = learnPMCode(m,p,nEpoch,nSave,respfile)
% LEARNPCODE learn pattern and motion bases (codes) in second layer
%
%   [m,p,snr] = learnPMCode(m,p,nEpoch,nSave,respfile)
%
% MooGu Z. <hzhu@case.edu>
% April 20, 2014 - Version 0.1

% number of iteration per save
if ~exist('nSave','var'), nSave = 10000; end

% Load First Layer Responds
if ~exist('respfile','var')
    [m,p,respfile] = cbaserespond(m,p);
end
load(respfile,'R1L');

% Initialize SNR
snr.logamp = zeros(nEpoch,1);
snr.dphase = zeros(nEpoch,1);

% Choose Test Set
tindex = randi(p.load_segments,130,1);

% Segment Index Function
segs = @(n) (n-1)*p.imszt+1;
sege = @(n) n*p.imszt;

fprintf('\nPattern-Bases Learning (%d Epoches) START @ %s\n',nEpoch,datestr(now));

for i = 1 : nEpoch
    I = randi(p.load_segments,nSave,1);
    for j = 1 : nSave
        % ==COLLECT DATA==
        if p.use_gpu
            P = gsingle(R1L.dPhase(:,segs(I(j))+1:sege(I(j))));
            A = gsingle(R1L.logAmp(:,segs(I(j)):sege(I(j))));
            valid = gsingle(1);
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
        % ==LEARN FORM CODE BASES==
        % Centralize and Scale
        A = bsxfun(@minus,A,m.loga_means);
        A = bsxfun(@times,A,m.loga_factors);
        % Infering and Adapting
        
        
    end
    % SNR Calculation
    
    % Save Model and Parameters
end

fprintf('Pattern-Bases Learning (%d Epoches) END @ %s\n',datestr(now));

end

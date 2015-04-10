% This script lauch transform-mask model learning process by automatically
% set corresponding optimization parameters according to current objective
% value

% FNAME TGTOBJ NPATTERN NTRANS SIGMA and MROOT MCODE are needed here

if ~exist('m','var')
    % initialize transform-mask model
    [m,rec] = tmmodel(video,'nepoch',30,'nadapt',7,'ninfer',7, ...
        'npattern',npattern,'ntrans',ntrans,'noiseprior',sigma.noise, ...
        'sparseprior',sigma.sparse,'slowprior',sigma.slow, ...
        'patternbasesmoothprior',sigma.smpat, ...
        'transformbasesmoothprior',sigma.smtrans, ...
        'NegativeCutProbability',0.3,'AlphaNormalization',false, ...
        'Verbose',2);
    
    % save model and video
    save([mroot,'/',mcode],'m','rec');
end

if ~exist('video','var')
    % load animation data
    if exist([fileparts(fname),'/tmdata.mat'],'file')
        load([fileparts(fname),'/tmdata.mat']);
    else
        [video.v,video.ffindex,video.res] = dataPrepare(fname);
    end
    % save animation
    save([mroot,'/',mcode],'video','-append');
end

% learning loop
while m.obj.value > tgtobj
    % optimization according to objective value
    if m.obj.value > 100 * tgtobj
        setANorm(m,false);
        m.ctrl.probNegCut = .5;
        [m,rec] = tmmodel(video,'model',m,'nepoch',30, ...
            'nadapt',13,'ninfer',13);
    elseif m.obj.value > 10 * tgtobj
        setANorm(m,false);
        m.ctrl.probNegCut = .7;
        [m,rec] = tmmodel(video,'model',m,'nepoch',17, ...
            'nadapt',30,'ninfer',30);
    elseif m.obj.value > 2 * tgtobj
        setANorm(m,false);
        m.ctrl.probNegCut = .9;
        [m,rec] = tmmodel(video,'model',m,'nepoch',13, ...
            'nadapt',70,'ninfer',70);
    else
        setANorm(m,true);
        m.ctrl.probNegCut = .95;
        [m,rec] = tmmodel(video,'model',m,'nepoch',7, ...
            'nadapt',70,'ninfer',70);
    end
    % save model and video
    save([mroot,'/',mcode],'m','video','rec');
end


function simulateMCode(C,I,model,param,varargin)
% ANIMATEMCODE produce animation of Motion Code learned from second level

% Input Parsering
p = inputParser;

p.addRequired('C',@(x) ismatrix(x)&&isnumeric(x));
p.addRequired('I',@(x) ismatrix(x)&&isnumeric(x));
p.addRequired('model',@isstruct);
p.addRequired('param',@isstruct);
p.addParamValue('Filename',datestr(now),@ischar);
p.addParamValue('frate',24,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('nframe',8,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('mfactor',.2,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));

p.parse(C,I,model,param,varargin{:});

fname   = ['./fig/',p.Results.Filename];
frate   = p.Results.frate;
nframe  = ceil(p.Results.nframe);
mfactor = p.Results.mfactor;

delta  = mfactor/nframe; % Delta of Motion Vector in each iteration

% Get first layer responds
Z = anim2z(repmat(I,1,2),model,param); Z = Z(:,1);
% Initialize Motion Code Animation
R = zeros(numel(Z),nframe*2);
R(:,1) = Z;

% Output Bases' animation one by one
[~,ncode]=size(C);
for c = 1 : ncode
    % Generate Motion Code Animation
    for i = 2 : nframe
        R(:,i) = R(:,i-1) .* exp(1j*delta*C(:,c));
    end
    R(:,nframe+1:nframe*2) = R(:,nframe:-1:1);
    anim = z2anim(R,model);
    % Write into GIF
    if ncode > 1
        gifname = [fname,'[',num2str(c),'].gif'];
    else
        gifname = [fname,'.gif'];
    end
    anim2gif(anim,[model.patch_sz,model.patch_sz],gifname,false,frate);
end

end

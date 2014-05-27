function animateBase(A,varargin)
% ANIMATEBASE produce animation of Base learned from first level

omodeSet = {'disp','gif'};
dmodeSet = {'comb','real'};

% Input Parsering
p = inputParser;

p.addRequired('A',@(x) ismatrix(x)&&isnumeric(x));
p.addOptional('drawmode','real',...
    @(x) any(validatestring(x,dmodeSet)));
p.addOptional('outmode','disp',...
    @(x) any(validatestring(x,omodeSet)));
p.addParamValue('Filename',datestr(now),@ischar);
p.addParamValue('T',2.67,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('nIteration',3,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('nframe',64,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('ncolor',256,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('resolution',256,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));

p.parse(A,varargin{:});

omode  = p.Results.outmode;
dmode  = p.Results.drawmode;
fname  = ['./fig/',p.Results.Filename];
t      = p.Results.T;
niter  = ceil(p.Results.nIteration);
nframe = ceil(p.Results.nframe);
ncolor = floor(p.Results.ncolor);
res    = floor(p.Results.resolution);

delay  = t/nframe;    % Time Delay for each frame
dPhi   = 2*pi/nframe; % Delta of Phi in each iteration

% Output Bases' animation one by one
[npixel,nbase]=size(A); sz=sqrt(npixel);
for b = 1 : nbase
    % Generate and Normalize Base Image
    Base = reshape(A(:,b),sz,sz);
    Base = Base / (2*max(abs(Base(:))));
    % Determine Scale Ratio
    [m,n] = size(Base);
    if (res^2 >= 4*m*n)
        scale = round(res/sqrt(m*n));
        Base  = kron(Base,ones(scale));
    end
    % Generate First Frame
    switch dmode
        case 'comb'
            [I,cmap] = rgb2ind(comp2img(Base),ncolor);
        case 'real'
            [I,cmap] = rgb2ind(repmat(real(Base)+0.5,[1,1,3]),ncolor);
        otherwise
            error('[AnimateBase] Undefined Draw Method!');
    end
    % Output
    switch lower(omode)
        case 'disp'
            f = figure(); axis image off;
            imshow(I,cmap); pause(delay);
            nframe = nframe * niter;
        case 'gif'
            if nbase == 1
                gifname = [fname,'.gif'];
            else
                gifname = [fname,'[',num2str(b),'].gif'];
            end
            imwrite(I,cmap,gifname,'gif','DelayTime',delay);
        otherwise
            error('[AnimateBase] Undefined Output Method!');
    end
    
    % Phase Changing
    for k = 1 : nframe-1
        % Generate Frame
        switch dmode
            case 'comb'
                I = rgb2ind(comp2img(Base*exp(-1j*dPhi*k)),cmap);
            case 'real'
                I = rgb2ind(repmat(real(Base*exp(-1j*dPhi*k))+0.5,[1,1,3]),...
                    cmap);
        end           
        % Output
        switch omode
            case 'disp'
                imshow(I,cmap); pause(delay);
            case 'gif'
                imwrite(I,cmap,gifname,'gif','WriteMode','append',...
                    'DelayTime',delay);
        end
    end
    
end

end
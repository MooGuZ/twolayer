function animateMCode(C,crds,varargin)
% ANIMATEMCODE produce animation of Motion Code learned from second level

omodeSet = {'disp','gif'};

% Input Parsering
p = inputParser;

p.addRequired('C',@(x) ismatrix(x)&&isnumeric(x));
p.addRequired('crds',@(x) ismatrix(x)&&isnumeric(x));
p.addOptional('outmode','disp',...
    @(x) any(validatestring(x,omodeSet)));
p.addParamValue('Filename',datestr(now),@ischar);
p.addParamValue('T',2.67,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('nframe',64,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('ncolor',256,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('resolution',512,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));
p.addParamValue('mfactor',.2,@(x) isscalar(x)&&isnumeric(x)&&isreal(x));

p.parse(C,crds,varargin{:});

omode   = p.Results.outmode;
fname   = ['./fig/',p.Results.Filename];
t       = p.Results.T;
nframe  = ceil(p.Results.nframe);
ncolor  = floor(p.Results.ncolor);
res     = floor(p.Results.resolution);
mfactor = p.Results.mfactor;

delay  = t/nframe;       % Time Delay for each frame
delta  = mfactor/nframe; % Delta of Motion Vector in each iteration
area   = [ ...           % Draw Area
    min(crds(1,:))-1,max(crds(1,:))+1, ...
    min(crds(2,:))-1,max(crds(2,:))+1];

% Initialize Figure Window
swGIF = false;
switch lower(omode)
    case 'disp'
        f = figure();
    case 'gif'
        f = figure('visible','off');
        swGIF = true;
end
set(f, ...
    'Color','k', ...
    'inverthardcopy','off', ...
    'position',[100,100,res,res]);

% Calculate Spacial Frequency of First Layer Bases
spfreq = sqrt(sum(crds(3:4,:).^2));
% Check Constant Component and Remove Corresponding Code Component
if any(spfreq==0)
    C(spfreq==0,:) = 0;
    spfreq = spfreq + eps;
end
% Motion Vector Direction for Each Bases
mvdir = bsxfun(@rdivide,crds(3:4,:),spfreq);

% Output Bases' animation one by one
[~,ncode]=size(C);
for c = 1 : ncode
    % Sorting and Normalization
    [~,index] = sort(abs(C(:,c)),'descend');
    % Motion Vectors
    mvamp = delta * C(index,c)' ./ spfreq(index);
    xmvec = mvamp .* mvdir(1,index);
    ymvec = mvamp .* mvdir(2,index);
    % Initial Position
    xpos = crds(1,index);
    ypos = crds(2,index);
    cval  = C(index,c)' / max(abs(C(:,c)));
    % Plot First Frame 
    scatter(xpos,ypos,30,cval,'filled'); 
    axis equal off; axis(area); drawnow
    if swGIF
        % Get image from plot
        [I,cmap] = rgb2ind(hardcopy(f,'-dzbuffer','-r0'),ncolor);
        % Create and Write into GIF File
        gifname = [fname,'[',num2str(c),'].gif'];
        imwrite(I,cmap,gifname,'gif','DelayTime',delay);
    end
    % Apply Motion to Each Bases
    for k = 1 : nframe-1
        xpos = xpos + xmvec;
        ypos = ypos + ymvec;
        % Plot current frame
        scatter(xpos,ypos,50,cval,'filled'); 
        axis equal off; axis(area); drawnow
        if swGIF
            I = rgb2ind(hardcopy(f,'-dzbuffer','-r0'),cmap);
            imwrite(I,cmap,gifname,'gif','WriteMode','append',...
                'DelayTime',delay);
        end
    end
    
end


end

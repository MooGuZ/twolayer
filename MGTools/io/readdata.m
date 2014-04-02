function anim = readdata(m,p,varargin)
% READDATA read data chunks according to setting
%  readdata(m,p) default read 1 patch from data directory according to M and P.

% Input Parsing
ps = inputParser;

ps.addRequired('m',@isstruct);
ps.addRequired('p',@isstruct);
ps.addParamValue('chunks',ceil(rand(1)*p.num_chunks),@isnumeric);
ps.addParamValue('frames',1:p.imszt,@isnumeric);
ps.addParamValue('patchsize',[m.patch_sz,m.patch_sz],@(x) isvector(x)&&isnumeric(x));
ps.addParamValue('position','random',@isvector);

ps.parse(m,p,varargin{:});

chunks = ps.Results.chunks;
frames = ps.Results.frames;
psize  = ps.Results.patchsize;

if ischar(ps.Results.position)
    vs = 1 + p.topmargin + p.BUFF;
    ve = 1 + p.imsz - p.BUFF - psize(1);
    hs = 1 + p.BUFF;
    he = 1 + p.imsz - p.BUFF - psize(2);
    pos = [floor(vs + rand(1)*(ve-vs)),floor(hs + rand(1)*(he-hs))];
else
    pos = ps.Results.position;
end

anim = zeros(prod(psize),numel(chunks)*numel(frames));
aind = 1; nframe = numel(frames);
for c = 1 : numel(chunks)
    % Read whole data chunk
    fid  = fopen([p.data_root,'/chunk',num2str(chunks(c))],'r','b');
    data = reshape(fread(fid,p.imsz*p.imsz*p.imszt,'float'),p.imsz,p.imsz,p.imszt);
    anim(:,aind:aind+nframe-1) = ...
        reshape(data(pos(1):pos(1)+psize(1)-1, pos(2):pos(2)+psize(2)-1, frames), ...
        prod(psize),nframe);
    aind = aind + nframe;
    fclose(fid);
end 

end
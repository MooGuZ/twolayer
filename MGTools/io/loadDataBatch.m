function Data = loadDataBatch(m,p)
% LOADDATABATCH would load data into frame vectors according to setting in
% paramters
%
%   Data = loadDataBatch(m,p)
%
% MooGu Z. <hzhu@case.edu>

% CHANGE LOG
% June 17, 2014 - Version 0.1
% Sept 23, 2015 - Version 0.2
%   1. add support to data of NPLab-3DMotion database
%   2. revise process of fvp to a more concise way, while replace RANDI with
%      - RANDPERM to generate samples more evenly distributed

switch p.data.type
    case 'fvp'
        Data = zeros(m.patch_sz^2,p.data.scope*p.data.nframe,'uint8');
        % generate a index list to indicate the order of input data to organized in
        % inner data structure. When data scope (number of data needed by program) is
        % greater than quantity of input data, ensure every input data loaded at least
        % once. To the opposite, pick a portion of input data to make each video loaded
        % at most once.
        if p.data.scope > p.data.quantity
            s = [randperm(p.data.quantity), randperm(p.data.quantity, p.data.scope - p.data.quantity)];
        else
            s = randperm(p.data.quantity,p.data.scope);
        end
        % Load GIF Files
        for i = 1 : p.data.scope
            Data(:,(i-1)*p.data.nframe+1:i*p.data.nframe) ...
                = im2uint8(gif2anim( ...
                    [p.data.path,p.data.nameList{s(i)}], ...
                    1:p.data.nframe));
        end
        
    case 'nplab-3dmotion'
        Data = zeros(m.patch_sz^2, p.data.scope * p.data.nframe, 'uint8');
        for i = 1 : p.data.scope
            % due to the structure that data stored under multiple layer of folders,
            % traverse is a quite complicate job to finish. Here, I use Monte-Carlo
            % method to fetch the sample randomly and as evenly as posible
            dpath = nplab3dRandPick(p.data.path);
            Data(:, (i-1) * p.data.nframe + 1 : i * p.data.nframe) ...
                = im2uint8(nplab3dLoadVideo(dpath, m, p));
        end
        
end

end
        

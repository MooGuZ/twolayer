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
%   3. combined 'fvp' and 'nplab-3dmotion' by force preprocessing of gathering
%      - video material in one folder to avoid overwhelm data load time

switch p.data.type
    case {'fvp', 'nplab-3dmotion'}
        Data = zeros(m.patch_sz^2,p.data.scope*p.data.nframe,'uint8');
        % generate a index list to indicate the order of input data to organized in
        % inner data structure. When data scope (number of data needed by program) is
        % greater than quantity of input data, ensure every input data loaded at least
        % once. To the opposite, pick a portion of input data to make each video loaded
        % at most once.
        findex = randperm(numel(p.data.nameList), p.data.scope);
        % load video one by one
        for i = 1 : p.data.scope
            Data(:, (i-1) * p.data.nframe + 1 : i * p.data.nframe) ...
                = im2uint8(gif2anim( ...
                    fullfile(p.data.path, p.data.nameList{findex(i)}), ...
                    1 : p.data.nframe));
        end
        
end

end
        

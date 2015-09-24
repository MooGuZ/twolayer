function N = nplab3dCount(p)
% NPLAB3DCOUNT return number of sequence stored under specified folder. 
% This program assume all the sequence are stored under 3 levels of
% folders as CONFIG - SCENE - SEQUENCE.
%
% MooGu Z. <hzhu@case.edu>

% CHANGE LOG
% ----------
% Sep 24, 2015 - Version 0.1 : initial commit with full functionality

N = 0;

% list CONFIG folders
cfgFolderList = subFolderList(p);
% dig into each CONFIG folder
for iConfig = 1 : numel(cfgFolderList)
    cfgPath = fullfile(p, cfgFolderList{iConfig});
    % list SCENE folders
    sceneFolderList = subFolderList(cfgPath);
    % dig into each SCENE folder
    for iScene = 1 : numel(sceneFolderList)
        scenePath = fullfile(cfgPath, sceneFolderList{iScene});
        % count each SEQUENCE folder
        N = N + numel(subFilderList(scenePath));
    end
end

end

% sub-routine to list sub-folder under specified path
function flist = subFolderList(p)
    % obtain file/folder mixed list
    flist = dir(p);
    % filtering to keep folders
    flist = {flist([flist(:).isdir]).name};
    % remove fake folder ('.' and '..')
    flist = flist(~ismember(flist, {'.', '..'}));
end
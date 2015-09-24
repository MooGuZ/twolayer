function seqPath = nplab3dRandPick(p)
% NPLAB3DRANDPICK randomly pick a folder containing video sequence data under specified path
%
% MooGu Z. <hzhu@case.edu>

% CHANGE LOG
% ----------
% Sept 24, 2015 - Version 0.1 : initial commit with full functionality

% list CONFIG folders
cfgFolderList = subFolderList(p);
% change to one of CONFIG folder randomly
cfgPath = fullfile(p, cfgFolderList{randi(numel(cfgFolderList))});
% list SCENE folders
sceneFolderList = subFolderList(cfgPath);
% enter one of SCENE folder randomly
scenePath = fullfile(cfgPath, sceneFolderList{randi(numel(sceneFolderList))});
% list SEQUENCE folders
seqFolderList = subFolderList(scenePath);
% choose one of SEQUENCE folder randomly
seqPath = fullfile(scenePath, seqFolderList{randi(numel(seqFolderList))});

end
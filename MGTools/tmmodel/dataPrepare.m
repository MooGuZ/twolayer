function [v,ffindex,res] = dataPrepare(fname)
% DATAPREPARE do read animation in specified folder or files to
% convert them to MATLAB data file format in the same folder for
% later processing. This process also contains operations that
% pre-process the data according to Transform-Mask model's
% requirement.
%
% MooGu Z. <hzhu@case.edu>
% Version 0.1 - Mar 5, 2015

% create file list according to type of input argument
switch exist(fname,'file')
  % single file
  case {2}
    flist = {fname};
  % folder
  case {7}
    flist = dir([fname,'/*.gif']);
    flist = strcat([fname,'/'],{flist(:).name});
  % unkown
  otherwise
    error('Input argument donot refer any folder or animation files!');
end

% read animation
ffindex = 1;
[v,res] = gif2anim(flist{1});
for i = 2 : numel(flist)
    ffindex = [ffindex,size(v,2)+1];
    v = [v,gif2anim(flist{i})];
end


% if this is a folder store video structure to file system
if exist(fname,'file') == 7
    video.ffindex = ffindex;
    video.res     = res;
    video.v       = v;
    save([fname,'/tmdata.mat'],'video');
end

end

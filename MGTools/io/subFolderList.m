function flist = subFolderList(p)
% list sub-folder under specified path

    % obtain file/folder mixed list
    flist = dir(p);
    % filtering to keep folders
    flist = {flist([flist(:).isdir]).name};
    % remove fake folder ('.' and '..')
    flist = flist(~ismember(flist, {'.', '..'}));
end
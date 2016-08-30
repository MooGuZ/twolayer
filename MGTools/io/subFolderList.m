function flist = subFolderList(p)
% list sub-folder under specified path
%
% Sept 28, 2015 - Updated function to ignore hidden folders

% obtain file/folder mixed list
flist = dir(p);
% create qualified sub-folder rules, ignore hidden and system folders
isSubFolder = @(f) f.isdir() && (f.name(1) ~= '.');
% find out all qualified sub-folders
flist = flist(arrayfun(isSubFolder, flist));

end
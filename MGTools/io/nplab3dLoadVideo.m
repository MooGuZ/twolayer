function frames = nplab3dLoadVideo(p, model, param)
% NPLAB3DLOADVIDEO read video data in images to inner data structure to contain video data of
% NPLab 3D Motion database.
%
% NOTE : this program would create a gif file as a cache version of image sequences.
%
% MooGu Z. <hzhu@case.edu>


[~, fname, ~] = fileparts(p);
cacheFileName = fullfile(p, [fname, '.gif']);
% when GIF cache file doesn't exist
if exist(cacheFileName, 'file') ~= 2
    imfile2gif(p, p, model.patch_sz, param.data.nframe);
end
% load GIF cache file into inner data structure
frames = gif2anim(cacheFileName);

end
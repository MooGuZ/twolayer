function imfile2gif(infd, outfd, frmsz, frmrt, orgfrmrt)
% This script would collect all the images file under specific folder to
% a gif file as an animation. And this function would recursively do the
% same collection work to the folders under the beginning folder.
%
% IMFILE2GIF(INFD, OUTFD) would collect all recursive collect all the image
% files under folder INFD to GIF files under folder OUTFD.
%
% IMFILE2GIF(INFD, OUTFD, FRAMESIZE, FRAMERATE) do the same thing as above 
% form. Besides would convert animate to specified frame size, FRAMESIZE,
% and scale time-axis to frame rate, FRAMERATE. We use 23.97 as the default
% frame rate to estimate the time duration of input image sequence. For
% example, if there are 36 images in the input folder, and specified output
% frame rate 12. Then the gif generated would have 18 frames. FRAMERATE is
% optional in this form.
%
% MooGu Z. <hzhu@case.edu>
% Apr 29, 2015

% CHANGE LOG
% ----------
% Sept 24, 2015 - updated CFUNC interface to only take one parameter
% Sept 25, 2015 - updated to assuming input video last for 1 second
%                 without specification.
% Sept 28, 2015 - added a switcher to suppress output information

if exist('frmsz','var')
    if exist('frmrt','var')
        if exist('orgfrmrt', 'var')
            cfunc = @(anim) vdsf(anim, frmsz, ...
                ceil(size(anim, 3) / orgfrmrt * frmrt));
        else
            cfunc = @(anim) vdsf(anim, frmsz, frmrt);
        end
    else
        cfunc = @(anim) vdsf(anim, frmsz, size(anim, 3));
    end
    collect(infd, outfd, '', cfunc);
else
    collect(infd, outfd);
end

end

function collect(infd, outfd, fname, cfunc)
% this is a helper function to collect image files under current folder
% and recursively collect images in folders under current folder.

showInfoVerbosely = false;
overwriteGif = false;
    

% image file extension name set
imExtSet = {'.jpg','.png','.bmp','.jpeg','.tiff'};

% default value
if ~exist('fname','var'), fname = ''; end

% check the contents in INFD
flist = dir(infd);
% initialize index for directory and image file
dirInd = false(1,numel(flist));
imfInd = false(1,numel(flist));
gifInd = false(1,numel(flist));
% filter the FLIST one by one
for i = 1 : numel(flist)
    % ignore hidden file and folders, including '.' and '..'
    if flist(i).name(1) == '.', continue; end
    % category each elements (folder/file) in current directory
    if flist(i).isdir
        dirInd(i) = true; % sub-folder
    else
       [~,~,ext] = fileparts(flist(i).name);
       if any(strcmpi(ext,imExtSet))
           imfInd(i) = true; % image file
       elseif strcmpi(ext, '.gif')
           gifInd(i) = true; % gif file
       end 
    end
end
% get name set of each type
dirSet = {flist(dirInd).name};
imfSet = {flist(imfInd).name};
gifSet = {flist(gifInd).name};

if ~isempty(gifSet) && ~overwriteGif
    % copy gif files to output folder in force (ignore exist file in output folder)
    cellfun(@(fname) copyfile(fname, outfd, 'f'), gifSet);
    
elseif ~isempty(imfSet)
    % show information
    if showInfoVerbosely
        fprintf('Working in Directory : %s \t', infd);
    end
        
    % fetch fname if necessary
    if isempty(fname)
        [~,fname,~] = fileparts(infd);
    end
    
    % get resolution information by first image
    im = imread(fullfile(infd,imfSet{1}));
    % check resolution
    imsize = size(im);
    % initialize matrix to store all images
    anim = zeros(imsize(1),imsize(2),numel(imfSet));
    % store first image
    if numel(imsize) > 2, im = rgb2gray(im); end
    anim(:,:,1) = im2double(im);
    % collect following image file into ANIM
    for i = 2 : numel(imfSet)
        im = imread(fullfile(infd,imfSet{i}));
        if numel(size(im)) > 2, im = rgb2gray(im); end
        anim(:,:,i) = im2double(im);
    end
    % convert animation to target format, if necessary
    if exist('cfunc','var')
        anim = cfunc(anim);
    end
    % check existence of output folder, create it if necessary
    if ~exist(outfd, 'file')
        mkdir(outfd);
    end
    % generate GIF file
    animsize = size(anim);
    anim2gif(reshape(anim,prod(animsize(1:2)),animsize(3)), ...
        fullfile(outfd, [fname,'.gif']));
    
    % show information
    if showInfoVerbosely
        disp('[DONE]');
    end
end

% append separator '-' for sub-folders, if necessary
if ~isempty(fname)
    fname = [fname,'-'];
end

% recursively collect images in sub-folders
if exist('cfunc','var')
    for i = 1 : numel(dirSet)
        collect(fullfile(infd,dirSet{i}),outfd, ...
            [fname,dirSet{i}],cfunc);
    end
else
    for i = 1 : numel(dirSet)
        collect(fullfile(infd,dirSet{i}),outfd, ...
            [fname,dirSet{i}]);
    end
end

end

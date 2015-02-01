% This script reform image frames of data into individual GIF files. And
% GIF files are named by seriers number of scene, focus point, and motion
% direction. In sum, there are 5 variables in this dataset.
%
%   D | Variable Discription
%   ------------------------
%   1 | Pixel Index
%   2 | Time Index
%   3 | Moving Direction Index
%   4 | Focus Point Index
%   5 | Scene Index
%
% Each movie frame have size MxN, then there are M*N pixels for each frame.
%
% This script should start at the root directory of images. And images are
% organized in directories in order of SCENE and FOCUSxSPEED.
%
% MooGu Z. <hzhu@case.edu>
% Jun 6, 2014

% Orginal Data Parameter
orgfps = 120;
tgtres = [32,32];
tgtfps = 24;
ndirct = 1;

% Write Message
note = ['This dataset is created by MooGu Z. @ ',datestr(now), ...
    ' with Elementary dataset created on Nov 25, 2014.'];

experimentPath = '/Users/MooGu/Dropbox/NPLab/201402-MotionPatternSeparation/data/Elementary/';
outputFilePath = [experimentPath,'FIG/downSampledByVDS/'];
outputFileNamePattern = [outputFilePath,'SCENE%dFOCUS%dDIR%d.gif'];

% Create output folder if necessary
if ~exist(outputFilePath,'dir')
	system(['mkdir -p "',outputFilePath,'"']);
end

% Read Image in Sequence
orgPwd = pwd;
cd(experimentPath);
scenelist = dir('./Scene *'); nscene = length(scenelist);
for i = 1 : nscene
    fprintf('Scene %d ',i);
    cd(scenelist(i).name);
    fdlist = dir('*F*D*');
    nfocus = length(fdlist) / ndirct;
    for j = 1 : nfocus
        for k = 1 : ndirct
            cd(fdlist((j-1)*ndirct+k).name);
            % Detect Pictures Format
            extname = '.jpg';
            frmlist = dir(['*',extname]);
            if isempty(frmlist)
                extname = '.bmp';
                frmlist = dir(['*',extname]);
            end
            nframe = length(frmlist);
            % Read First Frame
            img = imread(['0',extname]);
            img = rgb2ycbcr(img);
            img = im2double(img(:,:,1));
            orgres = size(img);
            MOV = zeros([orgres,nframe]);
            MOV(:,:,1) = img;
            % Collect Movie
            for f = 2 : nframe
                img = imread([num2str(f-1),extname]);
                img = rgb2ycbcr(img);
                img = im2double(img(:,:,1));
                MOV(:,:,f) = img;
            end
            % Down Sample by VDS
            nframe = round(nframe*(tgtfps/orgfps));
            MOV = vds(MOV,tgtres,nframe);
            % Save GIF Files
            fname = sprintf(outputFileNamePattern,i,j,k);
            anim2gif(reshape(MOV,prod(tgtres),nframe), ...
                tgtres,fname,false,tgtfps);
            fprintf('.');
            cd ..
        end
    end
    cd ..
    fprintf(' DONE\n');
end
cd(orgPwd);

% This scripts file run onebase analysis over a set of files and store
% corresponding results to well-organized files.

% Path Setting
videoPath = '/Users/MooGu/Dropbox/NPLab/201402-MotionPatternSeparation/data/Elementary/FIG/downSampledByVDSF';
figPath = '/Users/MooGu/Dropbox/NPLab/201402-MotionPatternSeparation/data/Elementary/FIG/onebaseAnalysis/VDSF';
% Operation Setting
nTestRound = 13;

% Create output folder if necessary
if ~exist(figPath,'dir')
    system(['mkdir -p "',figPath,'"']);
end
% Get file list
flist = dir([videoPath,'/*.gif']);
flist = {flist(:).name};
% Analyze Videos
for i = 1 : length(flist)
    % Get base name of video file
    fname = flist{i};
    fname = fname(1 : find(fname=='.',1)-1);
    
    bestTheta = 0;
    bestObjValue = Inf;
    for j = 1 : nTestRound
        disp([fname,' TestRound[',num2str(j),'] >>']);
        % Run onebase analysis
        [~,objValue,~,theta] = onebase([videoPath,'/',flist{i}]);
        % Comparing to current best one
        if objValue < bestObjValue
            bestObjValue = objValue;
            bestTheta = theta;
        end
    end
    % Run onebase optimiation initialized by best theta
    disp([fname,' Final Round >>' ]);
    [anim,objValue,base] = onebase([videoPath,'/',flist{i}],bestTheta);
    disp('===========================================================');
    
    % Save reconstructed animation
    anim2gif(anim,[32,32],[figPath,'/',fname,'-videoRec[', ...
        num2str(objValue),'].gif'],false);
    % Draw complex base plot
    dispBase(base,[figPath,'/',fname,'-basePlot']);
    % Save base animation
    animateBase(base,'outmode','gif','Filename', ...
        [figPath,'/',fname,'-baseAnim']);
end
    
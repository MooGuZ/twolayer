function [modelp,datap] = anaParamSetup(modelPath,dataPath,anaPath, ...
    batchSize,distBoxNum,unitOrder)
% ANAPARAMSETUP generate parameter structure of model and data preparing as the
% intput argument of function MODELANALYSIS.
%
%   [modelp,datap] = anaParamSetup(modelPath,dataPath,anaPath,batchSize,disBoxNum)
%   initialize parameter structure with corresponding information of model profiles
%   and dataset. In the arguments, MODELPATH and DATAPATH are necessary, others are
%   optional. If they are absent, default values would setup for them.
%
% see also modelAnalysis, anaPlot.
%
% MooGu Z. <hzhu@case.edu>
% -------------------------------------------------------------------------------------
% Version 0.1 [June 21, 2014] - Start Version
% Version 0.2 [June 23, 2014] - Add support for semilog distribution plots


modelp.path = modelPath;
if modelp.path(end) ~= '/'
    modelp.path = [modelp.path,'/'];
end
% search for model profile files
flist = dir([modelp.path,'*.mat']);
modelp.nameList = {flist(:).name};
% find number strings in model profile names
findNumStr = @(fname) regexp(fname,'-?\d+\.?\d*|-?\d*\.?\d+','match');
numStr = cellfun(findNumStr,modelp.nameList,'UniformOutput',false);
% calculate number of iterations for each model profile
niterCal = @(str) sum(cellfun(@str2num,str));
modelp.IDList = cellfun(niterCal,numStr);
% sort by number of iterations
[modelp.IDList,index] = sort(modelp.IDList,'ascend');
% reorganize nameList
modelp.nameList = modelp.nameList(index);

datap.path = dataPath;
if datap.path(end) ~= '/'
    datap.path = [datap.path,'/'];
end
% search for data files
flist = dir([datap.path,'*.gif']);
datap.nameList = {flist(:).name};

if exist('anaPath','var')
    if anaPath(end) == '/'
        datap.anaFilePath = anaPath;
    else
        datap.anaFilePath = [anaPath,'/'];
    end
else
    datap.anaFilePath = './';
end

if exist('batchSize','var')
    datap.batchSize = batchSize;
else
    datap.batchSize = 500;
end

if exist('distBoxNum','var')
    modelp.distBoxNum = distBoxNum;
else
    modelp.distBoxNum = 100;
end

if exist('unitOrder','var')
    modelp.unitOrder = unitOrder;
else
    modelp.unitOrder = -3;
end

load([modelp.path,modelp.nameList{1}]);
datap.npixel = m.patch_sz^2;
datap.nframe = p.data.nframe;
modelp.baseQuantity = m.N;

end

% fundamental setting
frmSize  = 32;
dataType = 'nplab-3dmotion';
dataPath = '/home/hxz244/data/NPLab3DMotionArchive';
% temperal functions
timestamp = @() upper(strrep(strrep(datestr(now), '-', ''), ' ', '-'));
% gpu module
p.use_gpu = gpuDeviceCount();
% initialization process
[m,p] = modelInitParameters(frmSize,dataType,dataPath);
if p.whitening.enable
    [m,p] = modelInitWhitening(m,p);
end
[m,p] = modelInitBases(m,p);
% overwrite parameters for test
p.data.scope = 10;
% autosave setting
p.autosave.path = fullfile('/home/hxz244/state/comodel/nplab3dtest', timestamp());
if ~isdir(p.autosave.path)
    mkdir(p.autosave.path);
end
% initialize random number generator
rng('shuffle');
% learning process
save_model([p.autosave.path,'init-',timestamp(),'.mat'],m,p);
[m,p] = learnCBases(m,p,3,10);
p.firstlayer.eta_dA_target = .25 * p.firstlayer.eta_dA_target;
[m,p] = learnCBases(m,p,3,10);
save_model([p.autosave.path,'final-',timestamp(),'.mat'],m,p);
% end
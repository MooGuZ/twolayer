function [m,p] = modelInitParameters(frmSize,dataType,dataPath)
% MODELINITPARAMETERS generate default running profile for Complex Bases Model in
% structure 'm' (model) and 'p' (parameter).
%
%   [m,p] = modelInitParamters(frmSize)
%
% MooGu Z. <hzhu@case.edu>
% June 17, 2014 - Version 0.1

m.patch_sz = frmSize;
% Estimate Number of Bases
m.N = frmSize^2;  % firstlayer basis functions
m.L = frmSize*5;  % phasetrans basis functions
m.K = frmSize*5;  % ampmodel basis functions

% specify priors
p.firstlayer.prior = 'slow_cauchy';
p.ampmodel.prior = 'slow_laplace';
p.phasetrans.prior = 'slow_cauchy';

% specify outerloop learning method
p.firstlayer.basis_method = 'steepest_adapt';
p.ampmodel.basis_method = 'steepest_adapt';
p.phasetrans.basis_method = 'steepest_adapt';

% specifiy inference methods
p.firstlayer.inference_method='minFunc_ind';
p.ampmodel.inference_method='minFunc_ind';
p.phasetrans.inference_method='minFunc_ind';

% data
p.data.type = dataType;
p.data.path = [dataPath,'/'];
if p.data.type = 'fvp'
    p.data.extname = '.gif';
    % Create Name List for Data Files
    flist = dir([p.data.path,'*',p.data.extname]);
    p.data.nameList = {flist(:).name}; clear flist
    p.data.quantity = numel(p.data.nameList);
    % The scale of data that processed at the same time
    p.data.scope = min(p.data.quantity,5000);
    p.data.nframe = 24;
end

% Save Path
p.autosave.path = './';

% Counter
m.t = zeros(3,1);

% Whitening
p.whitening.enable = true;
p.whitening.pixel_noise_fractional_variance = 0.01;
p.whitening.pixel_noise_variance_cutoff_ratio = 1.25;
p.whitening.X_noise_fraction = 8.0;
p.whitening.X_noise_var = 0.01;

% misc
p.use_gpu = false;      % this flag will use the GPU through the Jacket interface
p.renorm_length = 1;    % renorm basis function lengths after each iteration
p.normalize_crop = 0;   % normalize the data before inference
p.whiten_patches = 1;   % operate in the whitened domain
p.load_segments = 3200; % Responds collected from first layer
p.p_every = 0;          % display parameters
p.show_p = 0;           % display parameters
p.quiet = 0;            % display parameters

end
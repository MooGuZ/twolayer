function [Code,loglh,exitflag] = inferPCode(logAmp,m,p)
% INFERMCODE infer pattern code through the 2nd-layer of complex-bases model
%
% Usage: [Code,loglh,exitflag] = inferPCode(logAmp,m,p)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[~,N] = size(logAmp);

switch p.phasetrans.inference_method
  case {'minFunc_ind','minFunc_ind_lbfgs'}
    Code = .01 * m.B' * logAmp; % Pattern Code Initialization
    [Code,loglh,exitflag] = minFunc(@evalPCode,Code(:), ...
                                    p.ampmodel.minFunc_ind_Opts,logAmp,m,p);
    Code = reshape(Code,m.K,N);
end

end
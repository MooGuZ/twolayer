function [Code,loglh,exitflag] = inferMCode(dPhase,Mask,m,p)
% INFERMCODE infer motion code through the 2nd-layer of complex-bases model
%
% Usage: [Code,loglh,exitflag] = inferMCode(dPhase,Mask,m,p)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[~,N] = size(dPhase);

switch p.phasetrans.inference_method
  case {'minFunc_ind','minFunc_ind_lbfgs'}
    Code = .001 * m.D' * dPhase; % Motion Code Initialization
    [Code,loglh,exitflag] = minFunc(@evalMCode,Code(:), ...
                                    p.phasetrans.minFunc_ind_Opts,dPhase,Mask,m,p);
    Code = reshape(Code,m.L,N);
end

end
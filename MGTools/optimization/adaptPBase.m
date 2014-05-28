function [Base,noise,exitflag] = adaptPBase(Code,logAmp,m,p)
% INFERMCODE adapt pattern bases according by optimization methods.
%
% Usage: [Base,noise,exitflag] = adaptMBase(Code,dPhase,Mask,m,p)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[nbase,~] = size(Code);
    [N,~] = size(logAmp);

switch p.ampmodel.basis_method
  case {'minFunc_ind','minFunc_ind_lbfgs'}
    [Base,noise,exitflag] = minFunc(@evalPBase,m.B(:), ...
                                    p.ampmodel.minFunc_ind_Opts,Code,logAmp);
    m.B = reshape(Base,N,nbase);
end

end
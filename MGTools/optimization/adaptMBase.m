function [Base,loglh,exitflag] = adaptMBase(Code,dPhase,Mask,m,p)
% INFERMCODE adapt motion bases according by optimization method
%
% Usage: [Base,loglh,exitflag] = adaptMBase(Code,dPhase,Mask,m,p)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[nbase,~] = size(Code);
    [N,~] = size(dPhase);

switch p.phasetrans.basis_method
  case {'minFunc_ind','minFunc_ind_lbfgs'}
    [Base,loglh,exitflag] = minFunc(@evalMBase,m.D(:), ...
                                    p.phasetrans.minFunc_ind_Opts,Code,dPhase,Mask);
    m.D = reshape(Base,N,nbase);
end

end
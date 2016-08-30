function [m, p] = adaptPBase(Code,logAmp,m,p)
% INFERMCODE adapt pattern bases according by optimization methods.
%
% Usage: [Base,noise,exitflag] = adaptMBase(Code,dPhase,Mask,m,p)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[nbase,~] = size(Code);
    [N,~] = size(logAmp);

switch p.ampmodel.basis_method
  case {'steepset_adapt'}
    err = logAmp - m.B * Code;
    dB = err * Code';
    eta_dB = p.ampmodel.B_eta * dB / size(Code, 2);
    m.B = m.B + eta_dB;
    
    if max(abs(eta_dB(:))) > p.ampmodel.eta_dB_target
        p.ampmodel.B_eta = p.ampmodel.B_eta * p.ampmodel.down_factor;
    else
        p.ampmodel.B_eta = p.ampmodel.B_eta * p.ampmodel.up_factor;
    end
    
  case {'minFunc_ind','minFunc_ind_lbfgs'}
    Base = minFunc(@evalPBase,m.B(:), ...
                                    p.ampmodel.minFunc_ind_Opts,Code,logAmp);
    m.B = reshape(Base,N,nbase);
end

end
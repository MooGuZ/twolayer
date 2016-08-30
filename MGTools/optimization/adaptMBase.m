function [m, p] = adaptMBase(Code,dPhase,Mask,m,p)
% INFERMCODE adapt motion bases according by optimization method
%
% Usage: [Base,loglh,exitflag] = adaptMBase(Code,dPhase,Mask,m,p)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[nbase,~] = size(Code);
    [N,~] = size(dPhase);

switch p.phasetrans.basis_method
  case {'steepest_adapt'}
    err = Mask .* (dPhase - m.D * Code);
    dD = sin(err) * Code';
    eta_dD = p.phasetrans.D_eta * dD / size(Code, 2);
    m.D = m.D + eta_dD;
    
    if max(abs(eta_dD(:))) > p.phasetrans.eta_dD_target
        p.phasetrans.D_eta = p.phasetrans.D_eta * p.phasetrans.down_factor;
    else
        p.phasetrans.D_eta = p.phasetrans.D_eta * p.phasetrans.up_factor;
    end
    
  case {'minFunc_ind','minFunc_ind_lbfgs'}
    Base = minFunc(@evalMBase,m.D(:), ...
                                    p.phasetrans.minFunc_ind_Opts,Code,dPhase,Mask);
    m.D = reshape(Base,N,nbase);
end

end
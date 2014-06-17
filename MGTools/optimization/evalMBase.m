function [loglh,grad] = evalMBase(Base,Code,dPhase,Mask)
% EVALMCODE evaluate motion base by log-likelihood.
%
% Usage: [loglh,grad] = evalMCode(Base,Code,dPhase,Mask)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[nbase,~] = size(Code);
    [N,~] = size(dPhase);  
Base = reshape(Base,N,nbase);

delta = Mask .* (dPhase - Base * Code);

% Calculate Log-Likelihood
loglh = numel(delta)-sum(cos(delta(:))); % Noise(von Mises Distribution)

% Calculate Gradient for Log-Likelihood
if nargout > 1
    grad = -((Mask.*sin(delta))*Code');
    grad = grad(:);
end

end
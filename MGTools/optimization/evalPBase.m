function [noise,grad] = evalPBase(Base,Code,logAmp)
% EVALMCODE evaluate pattern base as noise energy.
%
% Usage: [noise,grad] = evalPBase(Base,Code,logAmp)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[nbase,~] = size(Code);
    [N,~] = size(logAmp);
Base = reshape(Base,N,nbase);

delta = logAmp - Base * Code;

% Calculate Noise Energy
noise = sum(delta(:).^2);

% Calculate Gradients of Log-Likelihood
if nargout > 1
    grad = -2*(delta*Code');
    grad = grad(:);
end

end

  
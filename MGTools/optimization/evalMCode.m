function [loglh,grad] = evalMCode(C,dPhase,Mask,m,p)
% EVALMCODE evaluate motion code as log-likelihood by noise energy, sparseness, and
% changing-slowness.
%
% Usage: [loglh,grad] = evalMCode(C,dPhase,Mask,m,p)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[~,N] = size(dPhase);  
C = reshape(C,m.L,N);

noise = Mask .* (dPhase - m.D*C);

% Calculate Log-Likelihood
dC = diff(C,1,2);
switch p.phasetrans.prior
  case 'slow_laplace'
    loglh = p.phasetrans.phase_noise_factor*(numel(noise)-sum(cos(noise(:)))) ... % Noise (von Mises)
            + p.phasetrans.w_laplace_beta*sum(abs(C(:))) ...                      % Sparse (Laplace)
            + p.phasetrans.w_lambda_S*.5*sum(sum(dC.^2));                         % Slow-Changing (Gaussian)
  case 'slow_cauchy'
    sigma = p.phasetrans.w_cauchy_sigma;
    loglh = p.phasetrans.phase_noise_factor*(numel(noise)-sum(cos(noise(:)))) ... % Noise (von Mises)
            + p.phasetrans.w_cauchy_beta*sum(log(1+(C(:)/sigma).^2)) ...          % Sparse (Cauchy)
            + p.phasetrans.w_lambda_S*.5*sum(sum(dC.^2));                         % Slow-Changing (Gaussian)
end

% Calculate Gradients for Log-Likelihood
if nargout > 1
    switch p.phasetrans.prior
      case 'slow_laplace'
        grad = -p.phasetrans.phase_noise_factor*(m.D'*(Mask.*sin(noise))) ...
               + p.phasetrans.w_laplace_beta*sign(C) ...
               + p.phasetrans.w_lambda_S*[-dC(:,1),-diff(dC,1,2),dC(:,end)];
      case 'slow_cauchy'
        grad = -p.phasetrans.phase_noise_factor*(m.D'*(Mask.*sin(noise))) ...
               + p.phasetrans.w_cauchy_beta*2*(C./(C.^2+sigma^2)) ...
               + p.phasetrans.w_lambda_S*[-dC(:,1),-diff(dC,1,2),dC(:,end)];
    end
    grad = grad(:);
end

end
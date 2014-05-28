function [loglh,grad] = evalPCode(C,logAmp,m,p)
% EVALMCODE evaluate pattern code as log-likelihood by noise energy, sparseness, and
% changing-slowness.
%
% Usage: [loglh,grad] = evalPCode(C,logAmp,m,p)
%
% MooGu Z. <hzhu@case.edu>
% May 27, 2014 - Version 0.1

[~,N] = size(logAmp);
C = reshape(C,m.K,N);

noise = logAmp - m.B * C;

% Calculate Log-Likelihood
dC = diff(C,1,2);
switch p.ampmodel.prior
  case 'slow_laplace'
    loglh = p.ampmodel.loga_noise_factor*.5*sum(noise(:).^2) ...       % Noise(Gaussian)
            + p.ampmodel.v_laplace_beta*sum(abs(C(:))) ...             % Sparse(Laplace)
            + p.ampmodel.v_lambda_S*.5*sum(dC(:).^2);                  % Slow-Changing(Gaussian)
  case 'slow_cauchy'
    sigma = p.ampmodel.v_cauchy_sigma;
    loglh = p.ampmodel.loga_noise_factor*.5*sum(noise(:).^2) ...       % Noise(Gaussian)
            + p.ampmodel.v_cauchy_beta*sum(log(1+(C(:)/sigma).^2)) ... % Sparse(Laplace)
            + p.ampmodel.v_lambda_S*.5*sum(dC(:).^2);                  % Slow-Changing(Gaussian)
end

% Calculate Gradients of Log-Likelihood
if nargout > 1
    switch p.ampmodel.prior
      case 'slow_laplace'
        grad = -p.ampmodel.loga_noise_factor*(m.B'*noise) ...
               + p.ampmodel.v_laplace_beta*sign(C) ...
               + p.ampmodel.v_lambda_S*[-dC(:,1),-diff(dC,1,2),dC(:,end)];
      case 'slow_cauchy'
        grad = -p.ampmodel.loga_noise_factor*(m.B'*noise) ...
               + p.ampmodel.v_cauchy_beta*2*(C./(C.^2+sigma^2)) ...
               + p.ampmodel.v_lambda_S*[-dC(:,1),-diff(dC,1,2),dC(:,end)];
    end
    grad = grad(:);
end

end

  
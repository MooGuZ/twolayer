function m = p2m(alpha,phi,beta,theta,bia)
% P2M reshape parameters in 4-D space to tighter structure in model

% get dimension info
npixel   = size(alpha,1);
npattern = size(alpha,2);
ntrans   = size(phi,3);
% reshape parameters
m.alpha  = reshape(alpha,[npixel,npattern]);
m.phi    = reshape(wrapToPi(phi),[npixel,ntrans]);
% reshape responds if exists
if exist('beta','var')
    nframe   = size(beta,4);
    m.beta   = reshape(beta,[npattern,ntrans,nframe]);
    m.theta  = reshape(wrapToPi(theta),[npattern,ntrans,nframe]);
    m.bia    = reshape(bia,[npattern,nframe]);
end

end
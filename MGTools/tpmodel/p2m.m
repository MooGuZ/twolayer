function m = p2m(alpha,phi,beta,theta,bia,m)
% P2M reshape parameters in 4-D space to tighter structure in model

% get dimension info
npixel   = size(alpha,1);
npattern = size(beta,2);
ntrans   = size(phi,3);
nframe   = size(theta,4);
% reshape parameters
m.alpha  = reshape(alpha,[npixel,npattern]);
m.phi    = reshape(phi,[npixel,ntrans]);
m.beta   = reshape(beta,[npattern,ntrans,nframe]);
m.theta  = reshape(theta,[npattern,ntrans,nframe]);
m.bia    = reshape(bia,[npattern,nframe]);

end
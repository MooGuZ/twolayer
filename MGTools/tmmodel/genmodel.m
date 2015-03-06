function v = genmodel(alpha,phi,beta,theta,bia)
% GENMODEL generate video from parameters of transform-mask model
%
% This implementation is hybrid with FOR LOOP and BSXFUN to achive the
% performance and follow memory limitation at the same time. According to
% the test, hybrid method get a better performance than pure BSXFUN and
% much better than FOR LOOP.
%
% MooGu Z. <hzhu@case.edu>
% Feb 18, 2015 - Version 0.2

npixel = size(alpha,1);
nframe = size(beta,4);
% initialize reconstructed video
v = zeros(npixel,1,1,nframe);
% reconstruct video frame by frame
for f = 1 : nframe
v(:,:,:,f) = sum(bsxfun(@times,alpha,bsxfun(@plus,bia(:,:,:,f), ...
    sum(bsxfun(@times,beta(:,:,:,f),cos(bsxfun(@minus,phi,theta(:,:,:,f)))),3))),2);
end

end

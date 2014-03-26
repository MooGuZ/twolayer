function Z = anim2z(frames,m,p)
% ANIM2Z convert animation in column frames into complex parameters infered
% by first layer

frames = m.whitenMatrix * bsxfun(@minus,frames,m.imageMean);

Z = infer_Z(frames,m,p);

end
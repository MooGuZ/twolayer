function anim = animateImage(patch,mcind,m,p)

nframe = 314;

z = anim2z(patch,m,p);
z = z(:,1);

phase = repmat(angle(z),1,nframe) + m.D(:,mcind) * (0:nframe-1);

znew = bsxfun(@times,exp(1j*phase),abs(z));

anim = z2anim(znew,m);

anim2gif(anim);

end
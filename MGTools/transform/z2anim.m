function frames = z2anim(Z,m)
% Z2I convert first layer parameters into animation in columnn frames

A = abs(Z); Phi = angle(Z);

frames = real(m.A)*(A.*cos(Phi)) + imag(m.A)*(A.*sin(Phi));

frames = bsxfun(@plus,m.dewhitenMatrix*frames,m.imageMean);

end
function [pcode,mcode] = z2pm(Z,m,p)
% Z2PC infers pattern and motion code from complex parameter matrix learned
% by first layer

warning('off','MATLAB:nearlySingularMatrix')

[~,nframe] = size(Z);
assert(nframe>1,'[Z2PM] Paramter Matrix Z should be an animation!');

a = abs(Z); phase = angle(Z);

pcode = infer_v(calc_logamp(a,m,p),m,p);

[dphase,avalind] = calc_dtphase(a,phase,m,p);
mcode = infer_w(dphase,avalind,m,p);

end
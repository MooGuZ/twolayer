function Z = pm2z(pcode,mcode,m)
% PM2Z compose pattern and motion code into complex parameters

a = exp(bsxfun(@plus,bsxfun(@rdivide,m.B*pcode,m.loga_factors),m.loga_means));

phase = m.D*mcode;
for i = 2 : size(phase,2)
    phase(:,i) = phase(:,i-1) + phase(:,i);
end

Z = a.*exp(1j*phase);

end
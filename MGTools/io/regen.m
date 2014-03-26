function anim = regen(patch,m,p)

z = anim2z(patch,m,p);

anim2gif(z2anim(z,m));

[pcode,mcode] = z2pm(z,m,p);

anim = z2anim(pm2z(pcode,mcode,m),m);

anim2gif(anim);

end
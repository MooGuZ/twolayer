function dispSecondlayer(m)
% DISPSECONDLAYER Display kernels of motion code and pattern code learning in second layer

% Display kernels of Pattern Codes
dispCode(m.B,m.Acoords,'PatternCode');
% Display kernels of Motion Codes
dispCode(m.D,m.Acoords,'MotionCode');

end
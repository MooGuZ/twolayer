% This file contains fragment scripts for debugging optimization part of
% transform-mask model

%% normalize alpha
normFactor = sum(abs(m.alpha),1) / m.ctrl.anorm;
m.alpha = bsxfun(@rdivide,m.alpha,normFactor);
m.bia   = bsxfun(@times,m.bia,normFactor);
m.beta  = bsxfun(@times,m.beta,normFactor);

%% check effectiveness of dAlpha
ss = 7e-2;
testAlpha = alpha - ss * dAlpha;
if ctrl.swNegCut
    testAlpha(testAlpha < 0) = 0;
    if ctrl.swANorm
        testAlpha = ...
            bsxfun(@rdivide,testAlpha,sum(abs(testAlpha),1)/ctrl.anorm);
    end
end
testDelta = v - genmodel(testAlpha,phi,beta,theta,bia);
testObj   = ...
    objFunc(testAlpha,phi,beta,theta,bia,testDelta, ...
    sigma,ffindex,resolution);
disp('=============================================================');
fprintf('Obj-Value   Improvement >> %9.2e\n',obj.value-testObj.value);
fprintf('Obj-Noise   Improvement >> %9.2e\n',obj.noise-testObj.noise);
fprintf('Obj-Sparse  Improvement >> %9.2e\n',obj.sparse-testObj.sparse);
fprintf('Obj-Slow    Improvement >> %9.2e\n',obj.slow-testObj.slow);
fprintf('Obj-SMPat   Improvement >> %9.2e\n',obj.smpat-testObj.smpat);
fprintf('Obj-SMTrans Improvement >> %9.2e\n',obj.smtrans-testObj.smtrans);

%% check effectiveness of dPhi
ss = 1e-5;
testPhi = alpha - ss * dPhi;
testDelta = v - genmodel(alpha,testPhi,beta,theta,bia);
testObj   = ...
    objFunc(alpha,testPhi,beta,theta,bia,testDelta, ...
    sigma,ffindex,resolution);
disp('=============================================================');
fprintf('Obj-Value   Improvement >> %5.2e\n',obj.value-testObj.value);
fprintf('Obj-Noise   Improvement >> %5.2e\n',obj.noise-testObj.noise);
fprintf('Obj-Sparse  Improvement >> %5.2e\n',obj.sparse-testObj.sparse);
fprintf('Obj-Slow    Improvement >> %5.2e\n',obj.slow-testObj.slow);
fprintf('Obj-SMPat   Improvement >> %5.2e\n',obj.smpat-testObj.smpat);
fprintf('Obj-SMTrans Improvement >> %5.2e\n',obj.smtrans-testObj.smtrans);

% END OF SCRIPT

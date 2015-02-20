function [beta,theta,bia,delta,obj,i] = ...
    inferGD(niter,alpha,phi,beta,theta,bia,delta,obj, ...
        v,sigma,ffindex,resolution,stepInit,stepMin)
if niter < 1 
    i = 0; return
end

step = stepInit;
for i = 1 : niter
    [dBeta,dTheta,dBia] = ...
        dCoefficient(alpha,phi,beta,theta,bia,delta,sigma,ffindex);
    newBeta  = beta - step * dBeta;
    newTheta = wrapToPi(theta - step * dTheta);
    newBia   = bia - step * dBia;
    newDelta = v - genmodel(alpha,phi,newBeta,newTheta,newBia);
    newObj   = ...
        objFunc(alpha,phi,newBeta,newTheta,newBia,newDelta, ...
            sigma,ffindex,resolution);
    while(newObj.value > obj.value)
        step = step / 2;
        if step < stepMin, break; end
        newBeta  = beta - step * dBeta;
        newTheta = wrapToPi(theta - step * dTheta);
        newBia   = bia - step * dBia;
        newDelta = v - genmodel(alpha,phi,newBeta,newTheta,newBia);
        newObj   = ...
            objFunc(alpha,phi,newBeta,newTheta,newBia,newDelta, ...
                sigma,ffindex,resolution);
    end
    if step < stepMin, break; end
    beta  = newBeta;
    theta = newTheta;
    bia   = newBia;
    delta = newDelta;
    obj   = newObj;
end
end

function [beta,theta,bia,delta,obj,i,ctrl] = ...
    inferGD(niter,alpha,phi,beta,theta,bia,delta,obj, ...
        sigma,ctrl,v,ffindex,res)
if niter < 1 
    i = 0; return
end

step = ctrl.inferInitStep;
for i = 1 : niter
    % calculate derivatives of coefficients
    [dBeta,dTheta,dBia] = ...
        dCoefficient(alpha,phi,beta,theta,bia,delta,sigma,ffindex);
    % find proper step size to gain improvement
    while (true)
        newBeta  = beta - step * dBeta;
        newTheta = theta - step * dTheta;
        newBia   = bia - step * dBia;
        newDelta = v - genmodel(alpha,phi,newBeta,newTheta,newBia);
        newObj   = ...
            objFunc(alpha,phi,newBeta,newTheta,newBia,newDelta, ...
            sigma,ffindex,res);
        % if no improvement, shrink the step size
        if newObj.value > obj.value
            step = step / 2;
            if step < ctrl.accuracy, return
            else continue
            end
        end
        % make the step
        beta  = newBeta;
        theta = newTheta;
        bia   = newBia;
        delta = newDelta;
        obj   = newObj;
        % move to next iteration
        break
    end
end
end

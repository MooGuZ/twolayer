function [alpha,phi,delta,obj,i] = ...
    adaptGD(niter,alpha,phi,beta,theta,bia,delta,obj, ...
        v,sigma,ffindex,resolution,stepInit,stepMin, ...
        swAlpha,swPhi)
if niter < 1 || ~(swAlpha || swPhi)
    i = 0; return
end

% randomlize switcher of negative cutter
swNegCut = rand(1) < 0.9;

step = stepInit;
for i = 1 : niter
    [dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia,delta, ...
        sigma,resolution,swAlpha,swPhi);
    if swAlpha
        newAlpha = alpha - step * dAlpha;
        % keep alpha positive
        if swNegCut, newAlpha(newAlpha < 0) = 0; end
    else
        newAlpha = alpha; 
    end
    if swPhi
        newPhi = wrapToPi(phi - step * dPhi);
    else
        newPhi = phi;
    end
    newDelta = v - genmodel(newAlpha,newPhi,beta,theta,bia);
    newObj   = ...
        objFunc(newAlpha,newPhi,beta,theta,bia,newDelta, ...
            sigma,ffindex,resolution);
    while(newObj.value > obj.value)
        step = step / 2;
        if step < stepMin, break; end
        if swAlpha
            newAlpha = alpha - step * dAlpha;
            % keep alpha positive
            if swNegCut, newAlpha(newAlpha < 0) = 0; end
        else
            newAlpha = alpha;
        end
        if swPhi
            newPhi = wrapToPi(phi - step * dPhi);
        else
            newPhi = phi;
        end
        newDelta = v - genmodel(newAlpha,newPhi,beta,theta,bia);
        newObj   = ...
            objFunc(newAlpha,newPhi,beta,theta,bia,newDelta, ...
                sigma,ffindex,resolution);
    end
    if step < stepMin, break; end
    alpha = newAlpha;
    phi   = newPhi;
    delta = newDelta;
    obj   = newObj;
end
end

function [alpha,phi,beta,theta,bia,delta,obj,i] = ...
    adaptGD(niter,alpha,phi,beta,theta,bia,delta,obj, ...
        v,sigma,ffindex,resolution,stepInit,stepMin, ...
        swAlpha,swPhi)

% deal with easy case
if niter < 1 || ~(swAlpha || swPhi)
    i = 0; return
end

% randomlize switcher of negative cutter
swNegCut = rand(1) < sigma.probNegCut;

step = stepInit;
for i = 1 : niter
    % calculate derivatives of bases
    [dAlpha,dPhi] = dBase(alpha,phi,beta,theta,bia,delta, ...
        sigma,resolution,swAlpha,swPhi);
    % step size search
    while(true)
        % update alpha
        if swAlpha
            newAlpha = alpha - step * dAlpha;
            % normalize alpha
            if sigma.swANorm
                newAlpha = bsxfun(@rdivide,newAlpha, ...
                    sum(abs(newAlpha),1)+eps) * sigma.anorm;
            end
            % keep alpha positive
            if swNegCut, newAlpha(newAlpha < 0) = 0; end
        else
            newAlpha = alpha;
        end
        % update phi
        if swPhi
            newPhi = phi - step * dPhi;
        else
            newPhi = phi;
        end
        % calculate new objective value
        newDelta = v - genmodel(newAlpha,newPhi,beta,theta,bia);
        newObj   = ...
            objFunc(newAlpha,newPhi,beta,theta,bia,newDelta, ...
                sigma,ffindex,resolution);
        % update step size when test step does not 
        % - improve objective value
        if newObj.value > obj.value
            step = step / 2;
            if step < stepMin, return
            else continue
            end
        end
        % make the step and scale parameters if negative-cut affective
        if sigma.swANorm && swNegCut
            normFactor = (sum(abs(newAlpha),1) + eps) / sigma.anorm;
            alpha = bsxfun(@rdivide,newAlpha,normFactor);
            phi   = newPhi;
            beta  = bsxfun(@times,beta,normFactor);
            bia   = bsxfun(@times,bia,normFactor);
            delta = newDelta;
            obj   = objFunc(alpha,phi,beta,theta,bia,delta, ...
                sigma,ffindex,resolution);
        else
            alpha = newAlpha;
            phi   = newPhi;
            delta = newDelta;
            obj   = newObj;
        end
        % go to next iteration
        break
    end
end
end
